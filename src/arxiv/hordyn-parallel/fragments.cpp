int main(int argc, char** argv) {
		/*
		 * Create MPI processes and get their process ids 
		 */
		MPI_Init(&argc, &argv);

		int num_processes, process_id = 0;

		MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
		MPI_Comm_rank(MPI_COMM_WORLD, &process_id);

		/*
		 * Calculate follicles . 
		 * If the number of follicles process is not equal, then abort
		 * program.
		 */
		const uint32_t num_follicles = 4;
		uint32_t num_follicles_per_process = 0;

		int error_flag = 0;
		if ((num_follicles % num_processes) != 0) {
				error_flag = 1;
		} 

		int error_global = 0;
		MPI_Allreduce(&error_flag, &error_global, 1, MPI_INT, MPI_LOR,
						MPI_COMM_WORLD);

		if (error_global) {
				if (process_id == 0) {
						std::cerr << "Number of follicles per process is not the same." << "\n";
				}
				MPI_Finalize();
				return 1;
		}

		num_follicles_per_process = num_follicles / num_processes;

		/*
		 *
		 */
		std::vector<uint32_t> follicle_ids;
		for (size_t index =0; index < num_follicles_per_process; index++) {
				follicle_ids.push_back(
								process_id * num_follicles_per_process + index
								);
		}

		/*
		 * Define output directory to write output files for follicle id
		 */
		// TODO: Remove serial 
		std::string directory_path = create_folder(
						"output/fullmodel/2012_gaussian/f", 0);

		std::vector<std::string> directory_path_vec;
		for (size_t index = 0; index < follicle_ids.size(); index++) {
				directory_path_vec.push_back(
								create_folder(
										"output/fullmodel/2012_gaussian/f",
									   	follicle_ids[index]
										)
								);
		}


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
		// TODO: Remove serial 
		BioParams bio_params;

		// TODO:  Declare follice biological parameters
		std::vector<BioParams> bio_params_vec;
		for (auto follicle_id : follicle_ids){
				if (follicle_id == 0)
						bio_params_vec.push_back(BioParams{});
				else if (follicle_id == 1) 
						bio_params_vec.push_back(BioParams{.tau_h=1.2});
				else 
						bio_params_vec.push_back(BioParams{.tau_h=0.7});
		}

		for (size_t fid = 0; fid < num_follicles_per_process; fid++) {
				std::cout << process_id << " " << bio_params_vec[fid].tau_h << "\n";
		}


		/*
		 * Define field variables to store follicle cell density,
		 * age velocity, maturation velocity, apoptosis, local and global
		 * maturity, local and global fsh
		 */
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
		std::vector<FollicleInitialDistribution> follicle_initial_distribution_vec;
		for (auto follicle_id : follicle_ids){
				if (follicle_id == 0)
						follicle_initial_distribution_vec.push_back(FollicleInitialDistribution{});
				else if (follicle_id == 1) 
						follicle_initial_distribution_vec.push_back(FollicleInitialDistribution{});
				else 
						follicle_initial_distribution_vec.push_back(FollicleInitialDistribution{});
		}

		compute_phi(domain, phi, follicle_initial_distribution_vec[0].cx,
					   	follicle_initial_distribution_vec[0].cy,
					   	follicle_initial_distribution_vec[0].sigma);


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


		MPI_Finalize();

		if (1 == 1) return 0;

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
				if (ONE_STEP) break;
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

