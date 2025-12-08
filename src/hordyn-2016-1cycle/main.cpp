#include "domain.hpp"
#include "utils.hpp"
#include "time.hpp"
#include "params.hpp"
#include "computations.hpp"
#include "compute_rhs.hpp"
#include <iostream>
#include <vector>
#include <mpi.h>


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int Np, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    const uint32_t Nf = 10;
    int base_follicles = Nf / Np; 
    int extra_follicles = Nf % Np;

    uint32_t my_Nf = (my_rank < extra_follicles)
        ?  base_follicles + 1
        : base_follicles;
    uint32_t my_start_index =  (my_rank < extra_follicles)
        ? my_rank * (base_follicles + 1)
        : extra_follicles * (base_follicles +1) + (my_rank - extra_follicles) * base_follicles;

    std::vector<uint32_t> my_follicle_ids;
    for (size_t i = 0; i < my_Nf; i++) {
        my_follicle_ids.push_back(my_start_index + i);
    }

    std::vector<std::string> dir_path(my_Nf);
    std::vector<Log> log(my_Nf);
    for (size_t i=0; i < my_Nf; i++) {
        dir_path[i] = create_folder("output/hordyn-2016/f", my_follicle_ids[i]);
        log[i] = Log {};
    }


    // Number of cycles
    //const uint32_t Nc = 1;
    // Duration of one cell cycle
    const double Dc = 1.0;
    const double ys = 0.3;
    const uint32_t cells_per_half_cycle = 60;
    const uint32_t Ny = 120;

    const Domain domain = Domain(Dc, ys, cells_per_half_cycle, Ny);

    std::vector<Time> time(my_Nf);
    for (size_t i = 0; i < my_Nf; i++) {
        time[i] = Time {};
        time[i].cfl = 0.4;
        time[i].end_time = 1.5;
    }

    std::vector<FollicleParams> follicle_params(my_Nf);
    for (size_t i = 0; i < my_Nf; i++) {
        follicle_params[i].gamma2 = 0.5 + my_follicle_ids[i] * 0.05;
    }

    for (size_t i = 0; i < my_Nf; i++) {
        std::cout << "My rank is " << my_rank
            << " follicle id is " << my_follicle_ids[i]
            << " gamma_2 is " << follicle_params[i].gamma2
            << "\n";
    }

    const uint32_t num_cols = domain.Nx;
    const uint32_t num_rows = domain.Ny;
    std::vector<std::vector<double>> phi(my_Nf, 
            std::vector<double>(num_rows * num_cols, 0.0));
    std::vector<std::vector<double>> g(my_Nf,
            std::vector<double>(num_rows * num_cols, 0.0));
    std::vector<std::vector<double>> h(my_Nf,
            std::vector<double>(num_rows * num_cols, 0.0));
    std::vector<std::vector<double>> source(my_Nf,
            std::vector<double>(num_rows * num_cols, 0.0));
    std::vector<double> mf(my_Nf, 0.0);
    std::vector<double> uf(my_Nf, 0.0);

    std::vector<double> phi_max(my_Nf, 0.0);
    std::vector<double> mass(my_Nf, 0.0);
    for (size_t i = 0; i < my_Nf; i++) {
        compute_phi(domain, phi[i], follicle_params[i]);
        phi_max[i] = find_abs_max_value(phi[i]);
        mass[i] = compute_mass(domain, phi[i]);
        mf[i] = compute_mf(domain, phi[i]);
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

    U = compute_U(global_bio_params, M);
    for (size_t i = 0; i < my_Nf; i++)  {
        uf[i] = compute_uf(follicle_params[i], mf[i], U);
        compute_g_h(domain, g[i], h[i], follicle_params[i], uf[i]);
        compute_source(domain, source[i], global_bio_params, U);
    }

    std::vector<double> dt(my_Nf);
    for (size_t i = 0; i < my_Nf; i++) {
        dt[i] = compute_dt(domain, time[i], g[i], h[i], source[i]);
    }

    double local_min_dt  = find_abs_min_value(dt);
    double global_min_dt = 0.0;
    MPI_Allreduce(&local_min_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    for (size_t i = 0; i < my_Nf; i++) {
        time[i].dt = global_min_dt;
    }

    for (size_t i = 0; i < my_Nf; i++) {
        write_matrix(dir_path[i] + "phi" + std::to_string(time[i].step_count) + ".txt", 
                phi[i], num_rows, num_cols);
        write_vector(dir_path[i]+"coords_x.txt", domain.coords_x);
        write_vector(dir_path[i]+"coords_y.txt", domain.coords_y);
        write_vector(dir_path[i]+"centers_x.txt", domain.centers_x);
        write_vector(dir_path[i]+"centers_y.txt", domain.centers_y);

        log[i].step_count.push_back(time[i].step_count);
        log[i].dt.push_back(time[i].dt);
        log[i].time.push_back(time[i].current_time);
        log[i].phi_max.push_back(phi_max[i]);
        log[i].mass.push_back(mass[i]);
        log[i].mf.push_back(mf[i]);
        log[i].M.push_back(M);
        log[i].U.push_back(U);
        log[i].uf.push_back(uf[i]);
    }

    // TODO: Remove priting to screen
    if (my_rank == 0) {
        for (size_t i = 0; i < my_Nf; i++) {
            std::cout
                << time[i].step_count
                << " dt " << time[i].dt  << " t " << time[i].current_time
                << " dt_g " << time[i].dt_g << " dt_h " << time[i].dt_h 
                << " dt_source " << time[i].dt_source << "\n"
                << " max_g " << time[i].max_g << " max_h " << time[i].max_h 
                << " max_l " << time[i].max_source << "\n"
                << " dx " << domain.dx << " dy " << domain.dy << "\n"
                << " phi_max " << phi_max[i] << " mass " << mass[i]
                << " m_f " << mf[i] << " M " << M 
                << " u_f " << uf[i] << " U " << U << "\n";
        }
    }
    std::vector<std::vector<double>> rhs1(my_Nf,
            std::vector<double>(num_rows * num_cols, 0.0));
    std::vector<std::vector<double>> rhs2(my_Nf,
            std::vector<double>(num_rows * num_cols, 0.0));
    std::vector<std::vector<double>> rhs3(my_Nf,
            std::vector<double>(num_rows * num_cols, 0.0));
    std::vector<std::vector<double>> phi_intermediate_1(my_Nf,
            std::vector<double>(num_rows * num_cols, 0.0));
    std::vector<std::vector<double>> phi_intermediate_2(my_Nf,
            std::vector<double>(num_rows * num_cols, 0.0));

    bool ONE_STEP = false;

    while (time[0].current_time < time[0].end_time) {
        for (size_t i = 0; i < my_Nf; i++) {
            rhs1[i] = compute_rhs(domain, time[i], phi[i],
                    g[i], h[i], source[i]);
            for (size_t j = 0; j < phi[i].size(); ++j) {
                phi_intermediate_1[i][j] = phi[i][j] - rhs1[i][j];
            }
            rhs2[i] = compute_rhs(domain, time[i], phi_intermediate_1[i],
                    g[i], h[i], source[i]);
            for (size_t j = 0; j < phi[i].size(); ++j) {
                phi_intermediate_2[i][j] = phi[i][j]
                    - (1.0/4.0) * rhs1[i][j]
                    - (1.0/4.0) * rhs2[i][j];
            }
            rhs3[i] = compute_rhs(domain, time[i], phi_intermediate_2[i],
                    g[i], h[i], source[i]);
            for (size_t j = 0; j < phi[i].size(); ++j) {
                phi[i][j] = phi[i][j]
                    - (1.0/6.0) * rhs1[i][j]
                    - (1.0/6.0) * rhs2[i][j]
                    - (2.0/3.0) * rhs3[i][j];
            }
            time[i].current_time  += time[i].dt;
            time[i].step_count++;

            phi_max[i] = find_abs_max_value(phi[i]);
            mass[i] = compute_mass(domain, phi[i]);
            mf[i] = compute_mf(domain, phi[i]);
        }

        local_mf_sum = 0.0;
        for (double val : mf) { local_mf_sum += val; }
        M = 0.0;
        MPI_Allreduce(&local_mf_sum, &M, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        U = 0.0;
        U = compute_U(global_bio_params, M);

        for (size_t i = 0; i < my_Nf; i++)  {
            uf[i] = compute_uf(follicle_params[i], mf[i], U);
            compute_g_h(domain, g[i], h[i], follicle_params[i], uf[i]);
            compute_source(domain, source[i], global_bio_params, U);
            dt[i] = compute_dt(domain, time[i], g[i], h[i], source[i]);
        }

        local_min_dt  = find_abs_min_value(dt);
        global_min_dt = 0.0;
        MPI_Allreduce(&local_min_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        for (size_t i = 0; i < my_Nf; i++) {
            time[i].dt = global_min_dt;
        }
        // Write at every time step
        for (size_t i = 0; i < my_Nf; i++) {
            write_matrix(dir_path[i] + "phi" + std::to_string(time[i].step_count) + ".txt", 
                    phi[i], num_rows, num_cols);
            write_vector(dir_path[i]+"coords_x.txt", domain.coords_x);
            write_vector(dir_path[i]+"coords_y.txt", domain.coords_y);
            write_vector(dir_path[i]+"centers_x.txt", domain.centers_x);
            write_vector(dir_path[i]+"centers_y.txt", domain.centers_y);

            log[i].step_count.push_back(time[i].step_count);
            log[i].dt.push_back(time[i].dt);
            log[i].time.push_back(time[i].current_time);
            log[i].phi_max.push_back(phi_max[i]);
            log[i].mass.push_back(mass[i]);
            log[i].mf.push_back(mf[i]);
            log[i].M.push_back(M);
            log[i].U.push_back(U);
            log[i].uf.push_back(uf[i]);
        }
        // TODO: Remove priting to screen
        if (my_rank == 0) {
            for (size_t i = 0; i < my_Nf; i++) {
                std::cout
                    << time[i].step_count
                    << " dt " << time[i].dt  << " t " << time[i].current_time
                    << " dt_g " << time[i].dt_g << " dt_h " << time[i].dt_h 
                    << " dt_source " << time[i].dt_source << "\n"
                    << " max_g " << time[i].max_g << " max_h " << time[i].max_h 
                    << " max_l " << time[i].max_source << "\n"
                    << " dx " << domain.dx << " dy " << domain.dy << "\n"
                    << " phi_max " << phi_max[i] << " mass " << mass[i]
                    << " m_f " << mf[i] << " M " << M 
                    << " u_f " << uf[i] << " U " << U << "\n";
            }
        }

        if (ONE_STEP) break;
    }

    for (size_t i = 0; i < my_Nf; i++) {
            log[i].write_summary(dir_path[i] + "summary.csv");
    }

    std::cout <<"End of program reached :), from rank "
        << my_rank <<  " follicles "
        << my_Nf <<"\n";

    MPI_Finalize();
    return 0;
}
