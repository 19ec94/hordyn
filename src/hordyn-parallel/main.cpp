#include "domain.hpp"
#include "params.hpp"
#include "computations.hpp"
#include "utils.hpp"
#include "compute_dt.hpp"
#include "compute_rhs.hpp"
#include "log.hpp"

#include <iostream>
#include <vector>
#include <stdexcept>
#include <mpi.h>


int main(int argc, char** argv) {
		MPI_Init(&argc, &argv);
		int num_processes, my_rank;
		MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

		const uint32_t total_num_follicles = 10;
		int base_follicles_per_process = total_num_follicles / num_processes;
		int extra_follicles = total_num_follicles % num_processes;

		uint32_t local_num_follicles = (my_rank < extra_follicles)
				? (base_follicles_per_process + 1)
				: base_follicles_per_process;

		uint32_t local_start_index = (my_rank < extra_follicles)
				? my_rank * (base_follicles_per_process + 1)
				: extra_follicles * (base_follicles_per_process + 1) + (my_rank - extra_follicles) * base_follicles_per_process;

		std::vector<uint32_t> local_follicle_ids;
		for (uint32_t i = 0; i < local_num_follicles; ++i) {
				local_follicle_ids.push_back(local_start_index + i);
		}


		std::vector<std::string> directory_path_vec(local_num_follicles);
		for (size_t i = 0; i < local_num_follicles; i++) {
				directory_path_vec[i] =
						create_folder( "output/fullmodel/parallel/10follicles/f",
										local_follicle_ids[i]);
		}

		std::vector<Log> log_vec(local_num_follicles);
		for (size_t i = 0; i < local_num_follicles; i++) {
				log_vec[i] = Log {};
		}

		//  Variables const for all follicles
		const uint32_t num_cycles = 8;
		const double single_cycle_duration = 1.0;
		const double threshold_y = 0.3;
		const uint32_t cells_per_half_cycle = 30;
		const uint32_t num_cells_y = 60;

		const Domain domain = Domain(num_cycles, single_cycle_duration,
						threshold_y, cells_per_half_cycle, num_cells_y);

		// Time parameters for each follicle
		std::vector<Time> time_vec(local_num_follicles);
		for (size_t  i = 0; i < local_num_follicles; i++){
				time_vec[i] = Time {};
				time_vec[i].cfl = 0.4;
				time_vec[i].end_time = 1.5;
		}

		std::vector<FollicleParams> follicle_params_vec(local_num_follicles);
		for (size_t i = 0; i < local_num_follicles; i++) {
				follicle_params_vec[i] = FollicleParams {};
				if (my_rank == 0) {
				follicle_params_vec[i].gamma2= 0.5 + i * 0.05;
				} else if (my_rank == 1) {
				follicle_params_vec[i].gamma2= 0.6 + i * 0.05;
				} else if (my_rank == 2) {
				follicle_params_vec[i].gamma2= 0.7 + i * 0.05;
				} else if (my_rank == 3) {
				follicle_params_vec[i].gamma2= 0.8 + i * 0.05;
				} else if (my_rank == 4) {
				follicle_params_vec[i].gamma2= 0.9 + i * 0.05;
				}
		}

		const uint32_t num_cols = domain.get_num_cols();
		const uint32_t num_rows = domain.get_num_rows();
		std::vector<std::vector<double>> phi(local_num_follicles, 
						std::vector<double>(num_rows * num_cols, 0.0));
		std::vector<std::vector<double>> g(local_num_follicles,
						std::vector<double>(num_rows * num_cols, 0.0));
		std::vector<std::vector<double>> h(local_num_follicles,
						std::vector<double>(num_rows * num_cols, 0.0));
		std::vector<std::vector<double>> source(local_num_follicles,
						std::vector<double>(num_rows * num_cols, 0.0));
		std::vector<double> mf(local_num_follicles, 0.0);
		std::vector<double> uf(local_num_follicles, 0.0);

		std::vector<double> phi_max_vec(local_num_follicles, 0.0);
		std::vector<double> mass_vec(local_num_follicles, 0.0);

		for (size_t i = 0; i < local_num_follicles; i++) {
				// Initialisation
				compute_phi(domain, phi[i],
								follicle_params_vec[i].cx,
								follicle_params_vec[i].cy,
								follicle_params_vec[i].sigma);

				// Computations
				phi_max_vec[i] = find_abs_max_value(phi[i]);
				mass_vec[i] = compute_mass(domain, phi[i]);
				mf[i] = compute_local_maturity(domain, phi[i]);
		}

		double local_mf_sum = 0.0;
		for (double val : mf) {
				local_mf_sum += val;
		}

		double M = 0.0;
		MPI_Allreduce(&local_mf_sum, &M, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		double U = 0.0;
		GlobalBioParams global_bio_params;
		global_bio_params.U_min = 0.9;
		global_bio_params.c = 10;
		global_bio_params.Lambda_bar = 1.0;

		U = compute_global_fsh(global_bio_params, M);

		for (size_t i = 0; i < local_num_follicles; i++)  {
				uf[i] = compute_local_fsh(follicle_params_vec[i], mf[i], U);
		}

		for (size_t i = 0; i < local_num_follicles; i++) {
				compute_g_h_source(domain, g[i], h[i], source[i],
								global_bio_params, follicle_params_vec[i], 
								uf[i], U);
		}

		std::vector<double> dt_vec(local_num_follicles);
		for (size_t i = 0; i < local_num_follicles; i++) {
				dt_vec[i] = compute_dt(domain, time_vec[i], g[i], h[i], source[i]);
		}

		double local_dt_min  = find_abs_min_value(dt_vec);
		double dt = 0.0;
		MPI_Allreduce(&local_dt_min, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		for (size_t i = 0; i < local_num_follicles; i++) {
				time_vec[i].dt = dt;
		}

		for (size_t i = 0; i < local_num_follicles; i++) {
				log_vec[i].write_matrix_to_file(
								directory_path_vec[i] + "phi" + std::to_string(time_vec[i].step_count) + ".txt", 
								phi[i], num_rows, num_cols);
				log_vec[i].write_vector_to_file(directory_path_vec[i]+"coords_x.txt", domain.coords_x);
				log_vec[i].write_vector_to_file(directory_path_vec[i]+"coords_y.txt", domain.coords_y);
				log_vec[i].write_vector_to_file(directory_path_vec[i]+"centers_x.txt", domain.centers_x);
				log_vec[i].write_vector_to_file(directory_path_vec[i]+"centers_y.txt", domain.centers_y);
				log_vec[i].step_count_vec.push_back(time_vec[i].step_count);
				log_vec[i].dt_vec.push_back(time_vec[i].dt);
				log_vec[i].time_vec.push_back(time_vec[i].current_time);
				log_vec[i].phi_max_vec.push_back(phi_max_vec[i]);
				log_vec[i].mass_vec.push_back(mass_vec[i]);
				log_vec[i].local_maturity_vec.push_back(mf[i]);
				log_vec[i].global_maturity_vec.push_back(M);
				log_vec[i].global_fsh_vec.push_back(U);
				log_vec[i].local_fsh_vec.push_back(uf[i]);
		}
		// TODO: Remove priting to screen
		if (my_rank == 0) {
				for (size_t i = 0; i < local_num_follicles; i++) {
						std::cout
								<< time_vec[i].step_count
								<< " dt " << time_vec[i].dt  << " t " << time_vec[i].current_time
								<< " dt_g " << time_vec[i].dt_g << " dt_h " << time_vec[i].dt_h 
								<< " dt_source " << time_vec[i].dt_source << "\n"
								<< " max_g " << time_vec[i].max_g << " max_h " << time_vec[i].max_h 
								<< " max_l " << time_vec[i].max_source << "\n"
								<< " dx " << domain.dx << " dy " << domain.dy << "\n"
								<< " phi_max " << phi_max_vec[i] << " mass " << mass_vec[i]
								<< " m_f " << mf[i] << " M " << M 
								<< " u_f " << uf[i] << " U " << U << "\n";
				}
		}

		std::vector<std::vector<double>> rhs1_vec (local_num_follicles,
						std::vector<double>(num_rows * num_cols, 0.0));
		std::vector<std::vector<double>> rhs2_vec (local_num_follicles,
						std::vector<double>(num_rows * num_cols, 0.0));
		std::vector<std::vector<double>> rhs3_vec (local_num_follicles,
						std::vector<double>(num_rows * num_cols, 0.0));
		std::vector<std::vector<double>> phi_star1_vec(local_num_follicles,
						std::vector<double>(num_rows * num_cols, 0.0));
		std::vector<std::vector<double>> phi_star2_vec(local_num_follicles,
						std::vector<double>(num_rows * num_cols, 0.0));

		while (time_vec[0].current_time < time_vec[0].end_time) {
				for (size_t i = 0; i < local_num_follicles; i++) {
						rhs1_vec[i] = compute_rhs(domain, time_vec[i], phi[i],
										g[i], h[i], source[i]);
						for (size_t j = 0; j < phi[i].size(); ++j) {
								phi_star1_vec[i][j] = phi[i][j] - rhs1_vec[i][j];
						}
						rhs2_vec[i] = compute_rhs(domain, time_vec[i], phi_star1_vec[i],
										g[i], h[i], source[i]);
						for (size_t j = 0; j < phi[i].size(); ++j) {
								phi_star2_vec[i][j] = phi[i][j]
										- (1.0/4.0) * rhs1_vec[i][j]
										- (1.0/4.0) * rhs2_vec[i][j];
						}
						rhs3_vec[i] = compute_rhs(domain, time_vec[i], phi_star2_vec[i],
										g[i], h[i], source[i]);
						for (size_t j = 0; j < phi[i].size(); ++j) {
								phi[i][j] = phi[i][j]
										- (1.0/6.0) * rhs1_vec[i][j]
										- (1.0/6.0) * rhs2_vec[i][j]
										- (2.0/3.0) * rhs3_vec[i][j];
						}
						time_vec[i].current_time  += time_vec[i].dt;
						time_vec[i].step_count++;

						phi_max_vec[i] = find_abs_max_value(phi[i]);
						mass_vec[i] = compute_mass(domain, phi[i]);
						mf[i] = compute_local_maturity(domain, phi[i]);
				}

				local_mf_sum = 0.0;
				for (double val : mf) {
						local_mf_sum += val;
				}
				M = 0.0;
				MPI_Allreduce(&local_mf_sum, &M, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				U = 0.0;
				U = compute_global_fsh(global_bio_params, M);

				for (size_t i = 0; i < local_num_follicles; i++)  {
						uf[i] = compute_local_fsh(follicle_params_vec[i], mf[i], U);
						compute_g_h_source(domain, g[i], h[i], source[i],
										global_bio_params, follicle_params_vec[i], 
										uf[i], U);
						dt_vec[i] = compute_dt(domain, time_vec[i], g[i], h[i], source[i]);
				}

				// new time step
				local_dt_min  = find_abs_min_value(dt_vec);
				MPI_Allreduce(&local_dt_min, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

				for (size_t i = 0; i < local_num_follicles; i++) {
						time_vec[i].dt = dt;
				}

				for (size_t i = 0; i < local_num_follicles; i++) {
						log_vec[i].write_matrix_to_file(
										directory_path_vec[i] + "phi" + std::to_string(time_vec[i].step_count) + ".txt", 
										phi[i], num_rows, num_cols);
						log_vec[i].write_vector_to_file(directory_path_vec[i]+"coords_x.txt", domain.coords_x);
						log_vec[i].write_vector_to_file(directory_path_vec[i]+"coords_y.txt", domain.coords_y);
						log_vec[i].write_vector_to_file(directory_path_vec[i]+"centers_x.txt", domain.centers_x);
						log_vec[i].write_vector_to_file(directory_path_vec[i]+"centers_y.txt", domain.centers_y);
						log_vec[i].step_count_vec.push_back(time_vec[i].step_count);
						log_vec[i].dt_vec.push_back(time_vec[i].dt);
						log_vec[i].time_vec.push_back(time_vec[i].current_time);
						log_vec[i].phi_max_vec.push_back(phi_max_vec[i]);
						log_vec[i].mass_vec.push_back(mass_vec[i]);
						log_vec[i].local_maturity_vec.push_back(mf[i]);
						log_vec[i].global_maturity_vec.push_back(M);
						log_vec[i].global_fsh_vec.push_back(U);
						log_vec[i].local_fsh_vec.push_back(uf[i]);
				}

				if (my_rank == 0) {
						for (size_t i = 0; i < local_num_follicles; i++) {
								std::cout
										<< time_vec[i].step_count
										<< " dt " << time_vec[i].dt  << " t " << time_vec[i].current_time
										<< " dt_g " << time_vec[i].dt_g << " dt_h " << time_vec[i].dt_h 
										<< " dt_source " << time_vec[i].dt_source << "\n"
										<< " max_g " << time_vec[i].max_g << " max_h " << time_vec[i].max_h 
										<< " max_l " << time_vec[i].max_source << "\n"
										<< " dx " << domain.dx << " dy " << domain.dy << "\n"
										<< " phi_max " << phi_max_vec[i] << " mass " << mass_vec[i]
										<< " m_f " << mf[i] << " M " << M 
										<< " u_f " << uf[i] << " U " << U << "\n";
						}
				}

				if (false) break;
		}
		for (size_t i = 0; i < local_num_follicles; i++) {
				log_vec[i].write_summary(directory_path_vec[i] + "summary.csv");
		}
		std::cout <<"End of program reached :), from rank " << my_rank <<  " follicles "<< local_num_follicles <<"\n";
		MPI_Finalize();
		return 0;
}

