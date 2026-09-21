#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

struct MeshParams {
    double m_min = 0.0;
    double m_max = 1.0;
    int Nm = 100000;
};


struct Grid {
    std::vector<double> centers;
    std::vector<double> edges;
    double dm;
    int Nm;
};

Grid setup_grid(const MeshParams& mesh_params) {
    Grid grid;
    grid.edges.resize(mesh_params.Nm + 1);
    grid.centers.resize(mesh_params.Nm);

    const double dm =
        (mesh_params.m_max - mesh_params.m_min) / mesh_params.Nm;

    grid.dm = dm;
    grid.Nm = mesh_params.Nm;

    for (int i = 0; i <= mesh_params.Nm; ++i) {
        grid.edges[i] = mesh_params.m_min + i * dm;
    }

    for (int i = 0; i < mesh_params.Nm; ++i) {
        grid.centers[i] = 0.5 * (grid.edges[i] + grid.edges[i + 1]);
    }

    return grid;
}

struct TimeParams {
    double T_end = 10.0;
    double CFL = 0.8;
    int output_frequency = 1;
};

struct BioParams {
    double xi = 15.297;
    double kappa_0 = 1.61786e-4;
    double delta = 0.0;
    double gamma_0 = 7.04742e-2;
    double eta = 5.58122e-3;
    double s_thr = 0.0;
    double nu = 2.47539;
    double alpha_0 = 0.033;
    double alpha_1 = 0.0;
    double alpha_2 = 0.0;
    double x_thr = 0.0;
};

struct Follicle {
    std::vector<double> phi_centers;
    double x;
    double s;
    double s_bar;
    std::vector<double> g_edges;
    double p_minus_lambda;

    Follicle(const MeshParams& mesh_params)
        : phi_centers(mesh_params.Nm, 0.0),
          x(0.0),
          s(0.0),
          s_bar(0.0),
          g_edges(mesh_params.Nm + 1, 0.0),
          p_minus_lambda(0.0) {}
};

struct InitCondition {
    struct MixtureComponent {
        double mu;
        double a;
        double b;
    };

    std::vector<MixtureComponent> components;

    InitCondition()
        : components({
              {27.5045, 0.0530822, 0.159247},
              {37.0671, 0.0376344, 0.112903},
              {49.7782, 0.0235043, 0.0705128},
          }) {}
};

void initialize_phi(
    const Grid& grid,
    const auto& init,
    std::vector<double>& phi_centers) {

    const double mu = init.mu;
    const double a = init.a;
    const double b = init.b;

    for (size_t j = 0; j < phi_centers.size(); ++j) {
        if (grid.centers[j] >= a && grid.centers[j] <= b) {
            phi_centers[j] = mu;
        } else {
            phi_centers[j] = 0.0;
        }
    }
}

void compute_macros(
    const Grid& grid,
    const std::vector<double>& phi_centers,
    double& x,
    double& s,
    double& s_bar) {

    x = 0.0;
    s = 0.0;

    for (size_t j = 0; j < phi_centers.size(); ++j) {
        x += phi_centers[j] * grid.dm;
        s += grid.centers[j] * phi_centers[j] * grid.dm;
    }

    s_bar = (x > 0.0) ? s / x : 0.0;
}


void compute_g(
    const Grid& grid,
    const BioParams& bio_params,
    const Follicle& follicle,
    std::vector<double>& g_edges) {

    const double alpha_0 = bio_params.alpha_0;
    const double xi = bio_params.xi;

    const double x = follicle.x;
    const double s_bar = follicle.s_bar;

    const double sigma = 0;//0.012;

    for (size_t i = 0; i < grid.edges.size(); ++i) {
        g_edges[i] =
            alpha_0 * (1.0 - grid.edges[i]) * x * s_bar
            - sigma * (xi - x) * (s_bar - grid.edges[i]);
    }
}

void compute_p_minus_lambda(
    const BioParams& bio_params,
    const std::vector<Follicle>& follicles,
    int i,
    double& p_minus_lambda) {

    const Follicle& follicle_i = follicles[i];

    const double x_i = follicle_i.x;
    const double xi = bio_params.xi;
    const double nu = bio_params.nu;
    const double gamma_0 = bio_params.gamma_0;
    const double kappa_0 = bio_params.kappa_0;
    const double eta_0 = bio_params.eta;

    const double prefactor = std::pow(xi - x_i, 1);

    double sum_other = 0.0;

    for (size_t j = 0; j < follicles.size(); ++j) {
        if (j != static_cast<size_t>(i)) {
            sum_other += std::pow(follicles[j].x, nu);
        }
    }

    const double term2 = gamma_0;
    const double term3_base = kappa_0;
    const double term3_inside =
        eta_0 * std::pow(x_i, nu) + sum_other;
    const double term3 = -term3_base * term3_inside;

    const double bracket = term2 + term3;

    p_minus_lambda = prefactor * bracket;
}

void compute_rhs(
    const Grid& grid,
    const BioParams& bio_params,
    std::vector<Follicle>& follicles,
    std::vector<std::vector<double>>& rhs_phi) {

    const int Nf = static_cast<int>(follicles.size());

    for (int i = 0; i < Nf; ++i) {
        compute_macros(
            grid,
            follicles[i].phi_centers,
            follicles[i].x,
            follicles[i].s,
            follicles[i].s_bar);

        compute_g(
            grid,
            bio_params,
            follicles[i],
            follicles[i].g_edges);

        compute_p_minus_lambda(
            bio_params,
            follicles,
            i,
            follicles[i].p_minus_lambda);
    }

    rhs_phi.assign(Nf, std::vector<double>(grid.Nm, 0.0));

    for (int i = 0; i < Nf; ++i) {
        auto& phi = follicles[i].phi_centers;
        auto& g = follicles[i].g_edges;
        auto& rhs = rhs_phi[i];

        const double pml = follicles[i].p_minus_lambda;

        // Boundary fluxes:
        // prescribed incoming PDF = 0,
        // outgoing flux uses the interior upwind value.

        // Left boundary: zero incoming PDF phi(0,t) = 0
        const double g_left = g[0];
        const double flux_left =
            (g_left < 0.0) ? g_left * phi[0] : 0.0;

        // Right boundary: zero incoming PDF phi(1,t) = 0
        const double g_right = g[grid.Nm];
        const double flux_right =
            (g_right > 0.0) ? g_right * phi[grid.Nm - 1] : 0.0;

        for (int j = 0; j < grid.Nm; ++j) {

            double F_left = flux_left;
            double F_right = flux_right;

            // Internal left interface
            if (j > 0) {
                const double gj = g[j];

                F_left =
                    (gj >= 0.0)
                        ? gj * phi[j - 1]
                        : gj * phi[j];
            }

            // Internal right interface
            if (j < grid.Nm - 1) {
                const double gj = g[j + 1];

                F_right =
                    (gj >= 0.0)
                        ? gj * phi[j]
                        : gj * phi[j + 1];
            }

            const double source = pml * phi[j];

            rhs[j] =
                (F_left - F_right) / grid.dm
                + source;
        }
    }
}


void rk4_step(
    const Grid& grid,
    const BioParams& bio_params,
    std::vector<Follicle>& follicles,
    double dt) {

    const int Nf = static_cast<int>(follicles.size());
    const int Nm = grid.Nm;

    std::vector<std::vector<double>> k1(
        Nf, std::vector<double>(Nm));
    std::vector<std::vector<double>> k2(
        Nf, std::vector<double>(Nm));
    std::vector<std::vector<double>> k3(
        Nf, std::vector<double>(Nm));
    std::vector<std::vector<double>> k4(
        Nf, std::vector<double>(Nm));
    std::vector<std::vector<double>> phi_backup(
        Nf, std::vector<double>(Nm));

    for (int i = 0; i < Nf; ++i) {
        phi_backup[i] = follicles[i].phi_centers;
    }

    compute_rhs(grid, bio_params, follicles, k1);

    for (int i = 0; i < Nf; ++i) {
        for (int j = 0; j < Nm; ++j) {
            follicles[i].phi_centers[j] =
                phi_backup[i][j] + 0.5 * dt * k1[i][j];
        }
    }

    compute_rhs(grid, bio_params, follicles, k2);

    for (int i = 0; i < Nf; ++i) {
        for (int j = 0; j < Nm; ++j) {
            follicles[i].phi_centers[j] =
                phi_backup[i][j] + 0.5 * dt * k2[i][j];
        }
    }

    compute_rhs(grid, bio_params, follicles, k3);

    for (int i = 0; i < Nf; ++i) {
        for (int j = 0; j < Nm; ++j) {
            follicles[i].phi_centers[j] =
                phi_backup[i][j] + dt * k3[i][j];
        }
    }

    compute_rhs(grid, bio_params, follicles, k4);

    for (int i = 0; i < Nf; ++i) {
        for (int j = 0; j < Nm; ++j) {
            follicles[i].phi_centers[j] =
                phi_backup[i][j]
                + dt * (
                    k1[i][j]
                    + 2.0 * k2[i][j]
                    + 2.0 * k3[i][j]
                    + k4[i][j]
                ) / 6.0;
        }
    }
}

int main() {
    const MeshParams mesh_params;
    const Grid grid = setup_grid(mesh_params);

    const BioParams bio_params;

    const int Nf = 3;
    std::vector<Follicle> follicles;
    follicles.reserve(Nf);

    for (int i = 0; i < Nf; ++i) {
        follicles.emplace_back(mesh_params);
    }

    InitCondition init_condition;

    for (int fol_idx = 0; fol_idx < Nf; ++fol_idx) {
        initialize_phi(
            grid,
            init_condition.components[fol_idx],
            follicles[fol_idx].phi_centers);
    }

    for (auto& follicle : follicles) {
        compute_macros(
            grid,
            follicle.phi_centers,
            follicle.x,
            follicle.s,
            follicle.s_bar);
    }

    for (auto& follicle : follicles) {
        compute_g(
            grid,
            bio_params,
            follicle,
            follicle.g_edges);
    }

    for (int i = 0; i < Nf; ++i) {
        compute_p_minus_lambda(
            bio_params,
            follicles,
            i,
            follicles[i].p_minus_lambda);
    }

    TimeParams time_params;
    double t = 0.0;

    std::filesystem::create_directories("cow-beta1-1");

    std::ofstream case_file("cow-beta1-1/case_log.csv");
    std::ofstream phi0_file("cow-beta1-1/phi0.csv");
    std::ofstream phi1_file("cow-beta1-1/phi1.csv");
    std::ofstream phi2_file("cow-beta1-1/phi2.csv");

    if (!case_file) {
        std::cerr << "Error: could not open case_log.csv for writing.\n";
        return 1;
    }

    if (!phi0_file) {
        std::cerr << "Error: could not open phi0.csv for writing.\n";
        return 1;
    }

    if (!phi1_file) {
        std::cerr << "Error: could not open phi1.csv for writing.\n";
        return 1;
    }

    if (!phi2_file) {
        std::cerr << "Error: could not open phi2.csv for writing.\n";
        return 1;
    }

    case_file << "t,dt,max_g,x0,s0,x1,s1,x2,s2\n";

    std::vector<std::vector<double>> phi0_history;
    std::vector<std::vector<double>> phi1_history;
    std::vector<std::vector<double>> phi2_history;

    phi0_history.push_back(follicles[0].phi_centers);
    phi1_history.push_back(follicles[1].phi_centers);
    phi2_history.push_back(follicles[2].phi_centers);

    int step = 0;

    while (t < time_params.T_end) {
        double max_g = 0.0;
        //double min_pml = 0.0;
        double max_pml = 0.0;

        for (const auto& follicle : follicles) {
            for (double gval : follicle.g_edges) {
                max_g = std::max(max_g, std::fabs(gval));
            }

            max_pml = std::max(max_pml,
                       std::fabs(follicle.p_minus_lambda));
        }

        double dt =
            time_params.CFL * grid.dm / (max_g + 1e-12);

        if (max_pml > 0.0) {
            const double dt_react = 0.5 / (max_pml + 1e-12);
            dt = std::min(dt, dt_react);
        }

        //dt = std::clamp(dt, 1e-6, 1e-2);

        dt = std::min(dt, time_params.T_end - t);
        rk4_step(grid, bio_params, follicles, dt);

        t += dt;
        ++step;

        if (step % time_params.output_frequency == 0) {
            std::cout
                << "t = " << std::setw(8) << t
                << "  dt = " << std::scientific << dt
                << " x0 = " << follicles[0].x
                << "  s0 = " << follicles[0].s
                << "  x1 = " << follicles[1].x
                << "  s1 = " << follicles[1].s
                << "  x2 = " << follicles[2].x
                << "  s2 = " << follicles[2].s
                << "\n";
        }

        if (step % time_params.output_frequency == 0) {
            case_file
                << t << ","
                << dt << ","
                << max_g << ","
                << follicles[0].x << ","
                << follicles[0].s << ","
                << follicles[1].x << ","
                << follicles[1].s << ","
                << follicles[2].x << ","
                << follicles[2].s
                << "\n";
        }

	if (step % time_params.output_frequency == 0) {
		std::cout
	    << "g(0) = " << follicles[0].g_edges[0]
	    << ", g(1) = " << follicles[0].g_edges[grid.Nm]
	    << '\n';
	}

        //if (step % time_params.output_frequency == 0) {
            //phi0_history.push_back(follicles[0].phi_centers);
            //phi1_history.push_back(follicles[1].phi_centers);
            //phi2_history.push_back(follicles[2].phi_centers);
        //}
    }

    auto write_phi_csv =
        [&](std::ofstream& file,
            const std::vector<std::vector<double>>& history) {

        file << "m";

        for (std::size_t k = 0; k < history.size(); ++k) {
            file << ",phi_t" << k;
        }

        file << "\n";

        for (int j = 0; j < grid.Nm; ++j) {
            file << grid.centers[j];

            for (const auto& snapshot : history) {
                file << "," << snapshot[j];
            }

            file << "\n";
        }
    };

    //write_phi_csv(phi0_file, phi0_history);
    //write_phi_csv(phi1_file, phi1_history);
    //write_phi_csv(phi2_file, phi2_history);

    case_file.close();
    phi0_file.close();
    phi1_file.close();
    phi2_file.close();

    return 0;
}
