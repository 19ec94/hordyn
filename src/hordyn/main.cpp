#include <iostream>
#include "domain.hpp"
#include "bio_params.hpp"
#include "computations.hpp"
#include "utils.hpp"
#include "compute_dt.hpp"
#include "compute_rhs.hpp"
#include "log.hpp"

#include <cmath>
#include <vector>


int main() {
		/*
		 * Define output directory to write output files
		 */
		int rank = 0;
		std::string directory_path = create_folder("output/fullmodel/2012_tau1.0_v2/f", rank);

		/*
		 * Instantiate the log variable to collect date over time steps
		 */
		Log log;

		/*
		 *  Create computational domain
		 */
		const uint32_t num_cycles = 8;
		const double single_cycle_duration = 1.0;
		const double threshold_y = 0.3;
		const uint32_t cells_per_half_cycle = 30;
		const uint32_t num_cells_y = 60;

		const Domain domain = Domain(num_cycles, single_cycle_duration,
						threshold_y, cells_per_half_cycle, num_cells_y);

		/*
		 * Define time parameters
		 */ 
		Time time;
		time.current_time = 0.0;
		time.cfl = 0.4;
		time.end_time = 8.0;
		

		/*
		 * Define biological parameters for follicles
		 */
		BioParams bio_params(threshold_y);

		/*
		 * Define field variables to store follicle cell density,
		 * age velocity, maturation velocity, apoptosis, local and global
		 * maturity, local and global fsh
		 */
		// const uint32_t num_follicles = 1;
		const uint32_t num_cols = domain.get_num_cols();
		const uint32_t num_rows = domain.get_num_rows();
		// Data are stored in row-major order
		std::vector<double> phi(num_rows * num_cols, 0.0);
		std::vector<double> g(num_rows * num_cols, 0.0);
		std::vector<double> h(num_rows * num_cols, 0.0);
		std::vector<double> source(num_rows * num_cols, 0.0);
		double local_fsh = 0.0;
		double local_maturity = 0.0;
		double global_fsh = 0.0;
		double global_maturity = 0.0;

		/*
		 * Initialize follicle cell density
		 */
		double cx = 0.30;
		double cy = 0.15;
		double sigma = std::sqrt(0.002); 
		compute_phi(domain, phi, cx, cy, sigma);

		/*
		 * Compute maximum follicular cell density, follicle mass, 
		 * local and global maturity, local and global fsh, 
		 * velocity functions g & h, and apoptosis
		 */
		double phi_max = 0.0;
		double mass = 0.0;
		phi_max = find_abs_max_value(phi);
		mass = compute_mass(domain, phi);
		local_maturity = compute_local_maturity(domain,phi);
		// TODO: MPI_SUM to calculate the global maturity
		global_maturity = local_maturity;
		global_fsh = compute_global_fsh(bio_params, global_maturity);
		local_fsh = compute_local_fsh(bio_params, local_maturity, global_fsh);
		compute_g_h_source(domain, g, h, source, bio_params, local_fsh, global_fsh);
		double dt = compute_dt(domain, time, g, h, source);
		// TODO: MPI_MIN to select the minimum dt of all the processes
		time.dt = dt;


		//std::vector<uint32_t> indices_doubling = domain.interface_x_index_doubling();
		//std::vector<uint32_t> indices_dirichlet = domain.cell_x_index_dirichlet();


		// TODO: Remove
		bool ONE_STEP = true;

		/* 
		 * Write the following data at time t=0:
		 * phi,
		 * domain.coords_x, domain.coords_y,
		 * domain.centers_x, domain.centers_y
		 * Log the following data at time t=0:
		 * step_count, dt, t,
		 * phi_max, mass
		 * local maturity, global maturity
		 * global fsh, local fsh
		 */
		log.write_matrix_to_file(
						directory_path + "phi" + std::to_string(time.step_count) + ".txt", 
						phi, num_rows, num_cols);
		log.write_vector_to_file(directory_path+"coords_x.txt", domain.coords_x);
		log.write_vector_to_file(directory_path+"coords_y.txt", domain.coords_y);
		log.write_vector_to_file(directory_path+"centers_x.txt", domain.centers_x);
		log.write_vector_to_file(directory_path+"centers_y.txt", domain.centers_y);
		log.step_count_vec.push_back(time.step_count);
		log.dt_vec.push_back(time.dt);
		log.time_vec.push_back(time.current_time);
		log.phi_max_vec.push_back(phi_max);
		log.mass_vec.push_back(mass);
		log.local_maturity_vec.push_back(local_maturity);
		log.global_maturity_vec.push_back(global_maturity);
		log.global_fsh_vec.push_back(global_fsh);
		log.local_fsh_vec.push_back(local_fsh);

		/*
		 * TODO: Print to screen. Remove it later.
		 */
		std::cout
				<< time.step_count
				<< " dt " << time.dt  << " t " << time.current_time
				<< " dt_g " << time.dt_g << " dt_h " << time.dt_h 
				<< " dt_source " << time.dt_source << "\n"
				<< " max_g " << time.max_g << " max_h " << time.max_h 
				<< " max_l " << time.max_source << "\n"
				<< " dx " << domain.dx << " dy " << domain.dy << "\n"
				<< " phi_max " << phi_max << " mass " << mass 
				<< " m_f " << local_maturity << " M " << global_maturity
				<< " u_f " << local_fsh << " U " << global_fsh << "\n";

		while (time.current_time < time.end_time ) {

				// Shu-Osher RK3
				std::vector<double> rhs1 (num_rows * num_cols, 0.0);
				std::vector<double> rhs2 (num_rows * num_cols, 0.0);
				std::vector<double> rhs3 (num_rows * num_cols, 0.0);
				std::vector<double> phi_star1(num_rows * num_cols, 0.0);
				std::vector<double> phi_star2(num_rows * num_cols, 0.0);

				rhs1 = compute_rhs(domain, time, phi, g, h, source);
				for (size_t i = 0;  i < phi.size(); ++i) {
						phi_star1[i] = phi[i] - rhs1[i];
				}
				rhs2 = compute_rhs(domain, time, phi_star1, g, h, source);
				for (size_t i = 0;  i < phi.size(); ++i) {
						phi_star2[i] = phi[i]
								- (1.0/4.0) * rhs1[i]
								- (1.0/4.0) * rhs2[i];
				}
				rhs3 = compute_rhs(domain, time, phi_star2, g, h, source);
				for (size_t i = 0;  i < phi.size(); ++i) {
						phi[i] = phi[i]
								- (1.0/6.0) * rhs1[i]
								- (1.0/6.0) * rhs2[i]
								- (2.0/3.0) * rhs3[i];
				}


				time.current_time += time.dt;
				time.step_count++;

				// Recalculate: phi_max, mass, local_maturity, global_maturity
				// global_fsh, local_fsh, g, h, source 
				phi_max = find_abs_max_value(phi);
				mass = compute_mass(domain, phi);
				local_maturity = compute_local_maturity(domain,phi);
				// TODO: MPI_SUM to gather all the follicles' local maturity
				global_maturity = local_maturity;
				global_fsh = compute_global_fsh(bio_params, global_maturity);
				local_fsh = compute_local_fsh(bio_params, local_maturity, global_fsh);
				compute_g_h_source(domain, g, h, source, bio_params, local_fsh, global_fsh);

				// recalculate time step
				dt = compute_dt(domain, time, g, h, source);
				// TODO: MPI_MIN to select the minimum of all processes
				time.dt = dt;

				/* 
				 * Write phi at time t
				 * Log the following data at time t: 
				 * step_count, dt, t,
				 * phi_max, mass
				 * local maturity, global maturity
				 * global fsh, local fsh
				 */
				log.write_matrix_to_file(
								directory_path + "phi" + std::to_string(time.step_count) + ".txt", 
								phi, num_rows, num_cols);
				log.step_count_vec.push_back(time.step_count);
				log.dt_vec.push_back(time.dt);
				log.time_vec.push_back(time.current_time);
				log.phi_max_vec.push_back(phi_max);
				log.mass_vec.push_back(mass);
				log.local_maturity_vec.push_back(local_maturity);
				log.global_maturity_vec.push_back(global_maturity);
				log.global_fsh_vec.push_back(global_fsh);
				log.local_fsh_vec.push_back(local_fsh);

				// TODO: Remove
				//if (time.current_time > 1.57) break;
				std::cout
						<< time.step_count
						<< " dt " << time.dt  << " t " << time.current_time
						<< " dt_g " << time.dt_g << " dt_h " << time.dt_h 
						<< " dt_source " << time.dt_source << "\n"
						<< " max_g " << time.max_g << " max_h " << time.max_h 
						<< " max_l " << time.max_source << "\n"
						<< " dx " << domain.dx << " dy " << domain.dy << "\n"
						<< " phi_max " << phi_max << " mass " << mass 
						<< " m_f " << local_maturity << " M " << global_maturity
						<< " u_f " << local_fsh << " U " << global_fsh << "\n";
		}
		/* 
		 * Write log variabels to a summary file.
		 */
		log.write_summary(directory_path + "summary.csv");
		std::cout << "End of program reached :)\n";
		return 0;
}
