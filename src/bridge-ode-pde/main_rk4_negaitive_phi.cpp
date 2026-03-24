#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <filesystem>
#include <algorithm>

struct MeshParams {
    double m_min = -0.1; // left boundary
    double m_max = +1.1; // right boundary
    int Nm = 5000; // Number of cells per follicle
};


struct Grid {
    std::vector<double> centers;
    std::vector<double> edges;
    double dm;
    int Nm;
};


Grid setup_grid(const MeshParams &mesh_params) {
    Grid grid;
    grid.edges.resize(mesh_params.Nm + 1);
    grid.centers.resize(mesh_params.Nm);

    const double dm = (mesh_params.m_max - mesh_params.m_min) / mesh_params.Nm;

    grid.dm = dm;
    grid.Nm = mesh_params.Nm;

    for (int i = 0; i < mesh_params.Nm + 1; ++i) {
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
    int output_frequency = 10;
};


struct BioParams {
    double xi = 16.0;
    double kappa_0 = 3.5;
    double delta = 0.1/16.0; //TODO: 16.0 -> xi
    double gamma_0 = 2.85;
    double eta = 4.22;
    double s_thr = 12.0;
    double nu = 2.0;
    double alpha_0 = 5.0;
    double alpha_1 = 2.0;
    double alpha_2 = 2.0;
    double x_thr = 12.0;
};


struct Follicle {
    // density
    std::vector<double> phi_centers;
    // mass 
    double x;
    // stage 
    double s;
    double s_bar;
    // maturation or stage velocity
    std::vector<double> g_edges;
    // cell doubling minus apoptosis rate
    double p_minus_lambda;

    Follicle (const MeshParams& mesh_params)
        : phi_centers(mesh_params.Nm, 0.0)
          , x(0.0)
          , s(0.0)
          , s_bar(0.0)
          , g_edges(mesh_params.Nm+1, 0.0)
          , p_minus_lambda(0.0)
    {}
};

struct InitCondition {
    struct MixtureComponent {
        double mu, a, b;
    };

    std::vector<MixtureComponent> components;

    InitCondition() : components({
            {50.0, 0.0100, 0.0900},
            {8.0,  0.0416, 0.2916}, 
            {60.0, 0.0130, 0.0764},
            {2.0,  0.8104, 0.8954}
            }) {}
};


void initialize_phi(const Grid& grid,
        const auto& init,
        std::vector<double>& phi_centers) {

    //const auto& comp = init_condition.components[fol_idx];
    double mu = init.mu;
    double a  = init.a;
    double b  = init.b;

    for (size_t j = 0; j < phi_centers.size(); ++j) {
        if (grid.centers[j] >= a && grid.centers[j] <= b) {
            phi_centers[j] = mu;
        } else {
            phi_centers[j] = 0.0;  // or whatever default you want
        }
    }

    return;
}


void compute_macros(const Grid& grid, const std::vector<double>& phi_centers,
        double &x, double &s, double &s_bar) {
    x = 0.0;
    s = 0.0;
    for (size_t j = 0; j < phi_centers.size(); ++j) {
        if (grid.centers[j] >= 0 && grid.centers[j] <= 1) {
            x += phi_centers[j] * grid.dm;
            s += grid.centers[j] * phi_centers[j] * grid.dm; 
        } 
    }
    s_bar = s / x;
}


void compute_g(const Grid& grid, const BioParams& bio_params,
        const Follicle& follicle,
        std::vector<double>& g_edges) {

    // inputs
    const double alpha_0 = bio_params.alpha_0;
    const double alpha_1 = bio_params.alpha_1;
    const double alpha_2 = bio_params.alpha_2;
    const double xi = bio_params.xi;
    const double x_thr = bio_params.x_thr;

    const double x = follicle.x;
    const double s_bar = follicle.s_bar;

    // Calculations
    double factor = std::pow( (x / xi) - s_bar, alpha_1);

    double hill_f = (
            std::pow(x, alpha_2))/
        (std::pow(x, alpha_2)+std::pow(x_thr, alpha_2)
        );

    double base = alpha_0 * factor * hill_f;
    for (size_t i = 0; i < grid.edges.size(); ++i) {
        g_edges[i] = base * (x - grid.edges[i] * xi);
    }
}


// compute source terms
void compute_p_minus_lambda(const BioParams& bio_params,
        const std::vector<Follicle>& follicles,
        int i,
        double& p_minus_lambda) {

    // Current follicle
    const Follicle& follicle_i = follicles[i];
    const double x_i = follicle_i.x;
    const double s_i = follicle_i.s;

    // Parameters
    const double xi = bio_params.xi;
    const double nu = bio_params.nu;
    const double delta = bio_params.delta;
    const double gamma_0 = bio_params.gamma_0;
    const double s_thr = bio_params.s_thr;
    const double kappa_0 = bio_params.kappa_0;
    const double eta_0 = bio_params.eta;  // assuming eta_0 = eta

    // Precompute common terms
    double prefactor = (1.0 / std::pow(xi, nu - 1.0)) * std::pow(xi - x_i, nu);

    // Sum over other follicles j≠i: x_j * s_j / ξ
    double sum_other = 0.0;
    for (size_t j = 0; j < follicles.size(); ++j) {
        if (j != static_cast<size_t>(i)) {
            sum_other += follicles[j].x * follicles[j].s / xi;
        }
    }

    // Bracketed term for each cell center
    //double m_k = grid.centers[k];  // cell center position

    // Individual terms in brackets
    double term1 = -delta;
    double term2 = gamma_0 * s_i / (s_thr + s_i);
    double term3_base = (kappa_0 / xi) * (1.0 - x_i / xi) * (x_i / xi);
    double term3_inside = eta_0 * (s_i / xi) * x_i + sum_other;
    double term3 = -term3_base * term3_inside;

    double bracket = term1 + term2 + term3;

    p_minus_lambda  = prefactor * bracket;
}


// compute RHS: dphi/dt = -(fluxR - fluxL)/dm + source
void compute_rhs(const Grid& grid,
                 const BioParams& bio_params,
                 std::vector<Follicle>& follicles,
                 std::vector<std::vector<double>>& rhs_phi) {
    int Nf = follicles.size();

    // update macros, g, p_minus_lambda for all follicles
    for (int i = 0; i < Nf; ++i) {
        compute_macros(grid, follicles[i].phi_centers, follicles[i].x, follicles[i].s, follicles[i].s_bar);
        compute_g(grid, bio_params, follicles[i], follicles[i].g_edges);
        compute_p_minus_lambda(bio_params, follicles, i, follicles[i].p_minus_lambda);
    }

    rhs_phi.assign(Nf, std::vector<double>(grid.Nm, 0.0));

    for (int i = 0; i < Nf; ++i) {
        auto& phi  = follicles[i].phi_centers;
        auto& g    = follicles[i].g_edges;
        auto& rhs  = rhs_phi[i];
        double pml = follicles[i].p_minus_lambda;

        for (int j = 0; j < grid.Nm; ++j) {
            double flux_left  = (j > 0) ? g[j]   * 0.5 * (phi[j-1] + phi[j])     : 0.0;
            double flux_right = (j < grid.Nm-1) ? g[j+1] * 0.5 * (phi[j]   + phi[j+1]) : 0.0;
            double source     = pml * phi[j];
            rhs[j] = ((-flux_right + flux_left) / grid.dm + source);
        }
    }
}

void rk4_step(const Grid& grid,
              const BioParams& bio_params,
              std::vector<Follicle>& follicles,
              double dt) {
    int Nf = follicles.size();
    int Nm = grid.Nm;

    // storage
    std::vector<std::vector<double>> k1(Nf, std::vector<double>(Nm));
    std::vector<std::vector<double>> k2(Nf, std::vector<double>(Nm));
    std::vector<std::vector<double>> k3(Nf, std::vector<double>(Nm));
    std::vector<std::vector<double>> k4(Nf, std::vector<double>(Nm));
    std::vector<std::vector<double>> phi_backup(Nf, std::vector<double>(Nm));

    // backup original phi
    for (int i = 0; i < Nf; ++i)
        phi_backup[i] = follicles[i].phi_centers;

    // k1
    compute_rhs(grid, bio_params, follicles, k1);

    // k2
    for (int i = 0; i < Nf; ++i)
        for (int j = 0; j < Nm; ++j)
            follicles[i].phi_centers[j] = phi_backup[i][j] + 0.5 * dt * k1[i][j];
    compute_rhs(grid, bio_params, follicles, k2);

    // k3
    for (int i = 0; i < Nf; ++i)
        for (int j = 0; j < Nm; ++j)
            follicles[i].phi_centers[j] = phi_backup[i][j] + 0.5 * dt * k2[i][j];
    compute_rhs(grid, bio_params, follicles, k3);

    // k4
    for (int i = 0; i < Nf; ++i)
        for (int j = 0; j < Nm; ++j)
            follicles[i].phi_centers[j] = phi_backup[i][j] + dt * k3[i][j];
    compute_rhs(grid, bio_params, follicles, k4);

    // combine
    for (int i = 0; i < Nf; ++i)
        for (int j = 0; j < Nm; ++j)
            follicles[i].phi_centers[j] =
                phi_backup[i][j]
                + dt * (k1[i][j] + 2.0*k2[i][j] + 2.0*k3[i][j] + k4[i][j]) / 6.0;
}



void time_step(const Grid& grid, const BioParams& bio_params, 
        std::vector<Follicle>& follicles, double dt) {

    int Nf = follicles.size();

    // 1. Update physics for ALL follicles
    for (int i = 0; i < Nf; ++i) {
        compute_macros(grid, follicles[i].phi_centers, follicles[i].x, follicles[i].s, follicles[i].s_bar);
        compute_g(grid, bio_params, follicles[i], follicles[i].g_edges);
        compute_p_minus_lambda(bio_params, follicles, i, follicles[i].p_minus_lambda);
    }

    std::vector<double> phi_new(grid.Nm);  // Copy current
                                           // 2. FVM update EACH follicle (Explicit Euler)
    for (int i = 0; i < Nf; ++i) {
        auto& phi = follicles[i].phi_centers;
        //std::vector<double> phi_new = phi;  // Copy current
        std::copy(phi.begin(), phi.end(), phi_new.begin());

        for (int j = 0; j < grid.Nm; ++j) {
            // Fluxes (0.5*(left+right))
            double flux_left = (j > 0) ? follicles[i].g_edges[j] * 0.5 * (phi[j-1] + phi[j]) : 0.0;
            double flux_right = (j < grid.Nm-1) ? follicles[i].g_edges[j+1] * 0.5 * (phi[j] + phi[j+1]) : 0.0;

            // Source
            double source = follicles[i].p_minus_lambda * phi[j];

            // FVM: φ_new = φ + dt*[-(F_R-F_L)/Δm + S]
            phi_new[j] = phi[j] + dt * ((-flux_right + flux_left) / grid.dm + source);
        }
        phi.swap(phi_new);  // Update
    }
}


int main() {
    // Define Mesh
    const MeshParams mesh_params;
    Grid grid = setup_grid(mesh_params);

    // Simulation parameters
    const BioParams bio_params;

    // Number of follicles
    int Nf = 4;
    std::vector<Follicle> follicles;
    follicles.reserve(Nf);
    for (int i = 0; i < Nf; ++i) {
        follicles.emplace_back(mesh_params);
    }

    InitCondition init_condition;

    for (int fol_idx =0; fol_idx < Nf; ++fol_idx) {
        initialize_phi(grid, init_condition.components[fol_idx], 
                follicles[fol_idx].phi_centers);
    }

    // Calculate mass (x) and stage (s)

    for (auto& follicle : follicles) {
        compute_macros(grid, follicle.phi_centers,
                follicle.x, follicle.s, follicle.s_bar);
        //std::cout << follicle.x <<  " " << follicle.s << " " << follicle.s_bar << "\n";
    }

    for (auto& follicle : follicles) {
        compute_g(grid, bio_params, follicle, follicle.g_edges);
    }


    for (int i = 0; i < Nf; ++i) {
        compute_p_minus_lambda(bio_params, follicles, 
                i, follicles[i].p_minus_lambda);
    }

    TimeParams time_params;
    double t = 0.0;

    std::ofstream case_file("case_log.csv");
    if (!case_file) {
        std::cerr << "Error: could not open case_log.csv for writing.\n";
        return 1;
    }
    std::ofstream phi0_file("phi0.csv");
    if (!phi0_file) {
        std::cerr << "Error: could not open phi0.csv for writing.\n";
        return 1;
    }
    std::ofstream phi1_file("phi1.csv");
    if (!phi1_file) {
        std::cerr << "Error: could not open phi1.csv for writing.\n";
        return 1;
    }
    std::ofstream phi2_file("phi2.csv");
    if (!phi2_file) {
        std::cerr << "Error: could not open phi2.csv for writing.\n";
        return 1;
    }
    std::ofstream phi3_file("phi3.csv");
    if (!phi3_file) {
        std::cerr << "Error: could not open phi3.csv for writing.\n";
        return 1;
    }


    case_file << "t,dt,max_g,x0,s0,x1,s1,x2,s2,x3,s3\n";

    static std::vector<std::vector<double>> phi0_history;  // store snapshots
    static std::vector<std::vector<double>> phi1_history;  // store snapshots
    static std::vector<std::vector<double>> phi2_history;  // store snapshots
    static std::vector<std::vector<double>> phi3_history;  // store snapshots
                                                           //
    phi0_history.push_back(follicles[0].phi_centers);  // assuming phi0 = follicle[0]
    phi1_history.push_back(follicles[1].phi_centers);  // assuming phi0 = follicle[0]
    phi2_history.push_back(follicles[2].phi_centers);  // assuming phi0 = follicle[0]
    phi3_history.push_back(follicles[3].phi_centers);  // assuming phi0 = follicle[0]
                                                           //
    int step = 0;
    while (t < time_params.T_end) {
        // recompute g first to estimate safe dt
        double max_g = 0.0;
        for (const auto& f : follicles)
            for (double gval : f.g_edges)
                max_g = std::max(max_g, std::fabs(gval));

        double dt = time_params.CFL * grid.dm / (max_g + 1e-12);
        dt = std::clamp(dt, 1e-6, 1e-2);


        //time_step(grid, bio_params, follicles, dt);
        rk4_step(grid, bio_params, follicles, dt);

        t += dt;
        step++;

        if (step % time_params.output_frequency == 0) {
        std::cout << "t = " << std::setw(8) << t 
            << "  dt = " << std::scientific << dt 
            << " x0 = " << follicles[0].x << "  "
            << " s0 = " << follicles[0].s << "  "
            << " x1 = " << follicles[1].x << "  "
            << " s1 = " << follicles[1].s << "  " 
            << " x2 = " << follicles[2].x << "  "
            << " s2 = " << follicles[2].s << "  "
            << " x3 = " << follicles[3].x << "  "
            << " s3 = " << follicles[3].s << "  "
            << "\n";
        }

        if (step % time_params.output_frequency == 0) {
            case_file << t << "," << dt << "," << max_g << ","
                << follicles[0].x << "," << follicles[0].s << ","
                << follicles[1].x << "," << follicles[1].s << ","
                << follicles[2].x << "," << follicles[2].s << ","
                << follicles[3].x << "," << follicles[3].s << "\n";
        }

        if (step % time_params.output_frequency == 0) {
            phi0_history.push_back(follicles[0].phi_centers);  // assuming phi0 = follicle[0]
            phi1_history.push_back(follicles[1].phi_centers);  // assuming phi0 = follicle[0]
            phi2_history.push_back(follicles[2].phi_centers);  // assuming phi0 = follicle[0]
            phi3_history.push_back(follicles[3].phi_centers);  // assuming phi0 = follicle[0]
        }

    }

    // helper lambda to write one phi history as CSV
    auto write_phi_csv = [&](std::ofstream& file,
            const std::vector<std::vector<double>>& history) {
        // header: m,phi_t0,phi_t1,...
        file << "m";
        for (std::size_t k = 0; k < history.size(); ++k)
            file << ",phi_t" << k;
        file << "\n";

        // data: one row per spatial point
        for (int j = 0; j < grid.Nm; ++j) {
            file << grid.centers[j];
            for (const auto& snapshot : history)
                file << "," << snapshot[j];
            file << "\n";
        }
    };

    write_phi_csv(phi0_file, phi0_history);
    write_phi_csv(phi1_file, phi1_history);
    write_phi_csv(phi2_file, phi2_history);
    write_phi_csv(phi3_file, phi3_history);

    case_file.close();
    phi0_file.close();
    phi1_file.close();
    phi2_file.close();
    phi3_file.close();

    case_file.close();

    return 0;
}

