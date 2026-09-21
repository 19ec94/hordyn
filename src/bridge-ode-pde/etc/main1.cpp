#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <filesystem>

class FvmSolver {
    public:
        struct Config {
            double m_min = -2.0; // Left boundary
            double m_max = 2.0; // Right boundary
            int Nm = 200;           // Number of cells per follicle
            int Nf = 3;             // Number of follicles
            double T = 2.0;         // Final time
            double CFL = 0.8;       // CFL number
            int output_freq = 20;   // Output every N steps
        };

        FvmSolver(const Config& cfg) : cfg_(cfg) {
            setup_grid();
            precompute_coefficients();
            initialize_solution();
            create_output_dir();
            write_initial_csv();
        }

        void solve_and_output() {
            print_header();

            double t = 0.0;
            int n = 0;
            while (t < cfg_.T - 1e-10) {
                double dt = compute_dt();
                std::cout << "dt " << dt << "\n";
                if (t + dt > cfg_.T) dt = cfg_.T - t;

                time_step(dt);
                t += dt;
                n++;

                if (n % cfg_.output_freq == 0 || t >= cfg_.T - 1e-10) {
                    write_csv(t, n);
                }

                if (n % 50 == 0) {
                    print_status(t, n);
                }
            }
            print_final_results();
        }

    private:
        Config cfg_;
        double dm_;
        std::vector<double> m_centers_, m_edges_, g_edges_, p_centers_, lambda_centers_;

        // Multi-follicle solution: phi[f][i] = density of follicle f at cell i
        std::vector<std::vector<double>> phi_, phi_new_;
        std::vector<std::vector<double>> F_;

        void setup_grid() {
            dm_ = (cfg_.m_max - cfg_.m_min) / cfg_.Nm;  // NEW: Δm = 4/Nm
            m_centers_.resize(cfg_.Nm);
            m_edges_.resize(cfg_.Nm + 1);
            g_edges_.resize(cfg_.Nm + 1);

            // NEW: Symmetric grid [-2, 2]
            for (int i = 0; i <= cfg_.Nm; ++i) {
                m_edges_[i] = cfg_.m_min + i * dm_;
            }
            for (int i = 0; i < cfg_.Nm; ++i) {
                m_centers_[i] = 0.5 * (m_edges_[i] + m_edges_[i + 1]);
            }
        } 

        void precompute_coefficients() {
            p_centers_.resize(cfg_.Nm);
            lambda_centers_.resize(cfg_.Nm);
            double center = 0.0;

            for (int i = 0; i < cfg_.Nm; ++i) {
                double m = m_centers_[i];
                double gaussian_8 = std::exp(-8.0 * (m - center) * (m - center));
                lambda_centers_[i] = 0.8 * gaussian_8;
                p_centers_[i] = 0.2 * gaussian_8;
            }

            for (int i = 0; i <= cfg_.Nm; ++i) {
                double m = m_edges_[i];
                g_edges_[i] = 0.5 * std::exp(-10.0 * (m - center) * (m - center));
            }
        }

        void initialize_solution() {
            phi_.resize(cfg_.Nf, std::vector<double>(cfg_.Nm));
            phi_new_.resize(cfg_.Nf, std::vector<double>(cfg_.Nm));
            F_.resize(cfg_.Nf, std::vector<double>(cfg_.Nm + 1));
            // NEW: Exponential decay parameters per follicle
            double mu[] = {2.0, 1.5, 1.0};      // Peak densities
            double sigma[] = {5.0, 4.0, 3.0};   // Decay rates

            // Characteristic window [0.1, 0.2]
            const double m_low = 0.1;
            const double m_high = 0.2;

            for (int f = 0; f < cfg_.Nf; ++f) {
                for (int i = 0; i < cfg_.Nm; ++i) {
                    double m = m_centers_[i];
                    // Sharp version
                    // Characteristic function: 1 if m ∈ [0.1, 0.2], else 0
                    //double chi = (m >= m_low && m <= m_high) ? 1.0 : 0.0;

                    //phi_[f][i] = mu[f] * std::exp(-sigma[f] * m) * chi;

                    // Smooth version
                    double chi_smooth = 0.5 * (1.0 + std::tanh(50.0 * (m - m_low))) * 0.5 * (1.0 - std::tanh(50.0 * (m - m_high)));
                    phi_[f][i] = mu[f] * std::exp(-sigma[f] * m) * chi_smooth;
                    //if (m >= 0.0 && m <= 1.0) {  // Restrict to [0,1] biology region
                    //    phi_[f][i] = mu[f] * std::exp(-sigma[f] * m);
                    //} else {
                    //    phi_[f][i] = 0.0;  // Zero outside [0,1]
                    //}
                }
            }
        }

        double compute_dt() const {
            double g_max = 0.0;
            for (double g : g_edges_) {
                g_max = std::max(g_max, std::abs(g));
            }
            return cfg_.CFL * dm_ / (g_max > 1e-12 ? g_max : 1.0);
        }

        void compute_fluxes() {
            for (int f = 0; f < cfg_.Nf; ++f) {
                F_[f][0] = 0.0;  // No inflow

                for (int i = 0; i < cfg_.Nm; ++i) {
                    if (g_edges_[i + 1] >= 0.0) {
                        F_[f][i + 1] = g_edges_[i + 1] * phi_[f][i];
                    } else {
                        F_[f][i + 1] = g_edges_[i + 1] * phi_[f][i + 1];
                    }
                }

                //F_[f][cfg_.Nm] = g_edges_[cfg_.Nm] * phi_[f][cfg_.Nm - 1];
                F_[f][cfg_.Nm] = 0.0; // No mass leaves
            }
        }

        void time_step(double dt) {
            compute_fluxes();
            double inv_dm_dt = dt / dm_;

            for (int f = 0; f < cfg_.Nf; ++f) {
                for (int i = 0; i < cfg_.Nm; ++i) {
                    double flux_div = (F_[f][i + 1] - F_[f][i]) * inv_dm_dt;
                    double source = dt * (p_centers_[i] - lambda_centers_[i]) * phi_[f][i];
                    phi_new_[f][i] = phi_[f][i] - flux_div + source;
                    phi_new_[f][i] = std::max(0.0, phi_new_[f][i]);
                }
            }
            phi_ = phi_new_;
        }

        // ✅ CORRECTED: Individual mass per follicle
        std::vector<double> masses() const {
            std::vector<double> masses(cfg_.Nf);
            for (int f = 0; f < cfg_.Nf; ++f) {
                masses[f] = std::reduce(phi_[f].begin(), phi_[f].end(), 0.0) * dm_;
            }
            return masses;
        }

        double total_mass_all() const {
            auto ms = masses();
            return std::reduce(ms.begin(), ms.end(), 0.0);
        }

        void create_output_dir() {
            std::filesystem::create_directory("output");
        }

        void write_csv(double t, int step) {
            std::ofstream file("output/phi_t" + std::to_string(step) + ".csv");
            file << "m";
            for (int f = 0; f < cfg_.Nf; ++f) {
                file << ",phi_f" << f;
            }
            file << ",g,lambda,p\n";

            for (int i = 0; i < cfg_.Nm; ++i) {
                double m = m_centers_[i];
                file << m;
                for (int f = 0; f < cfg_.Nf; ++f) {
                    file << "," << phi_[f][i];
                }
                double g_val = 0.5 * std::exp(-10.0 * (m - 0.3) * (m - 0.3));
                file << "," << g_val << "," << lambda_centers_[i] << "," << p_centers_[i] << "\n";
            }
            file.close();
        }

        void write_initial_csv() {
            write_csv(0.0, 0);
        }

        void print_header() const {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "\n=== Multi-Follicle 1D FVM Model ===\n";
            std::cout << "Nm = " << cfg_.Nm << ", Nf = " << cfg_.Nf 
                << ", T = " << cfg_.T << ", CFL = " << cfg_.CFL << "\n";
            std::cout << "CSV: m,phi_f0,phi_f1,phi_f2,g,lambda,p\n\n";
            std::cout << "Step\tTime\tMass_f0\t\tMass_f1\t\tMass_f2\t\tTotal\n";
            std::cout << "------------------------------------------------\n";
        }

        // ✅ CORRECTED: Print individual masses
        void print_status(double t, int n) const {
            auto ms = masses();
            std::cout << std::setw(4) << n << "\t" << std::setw(7) << t;
            for (int f = 0; f < cfg_.Nf; ++f) {
                std::cout << "\t" << std::setw(9) << ms[f];
            }
            std::cout << "\t" << std::setw(9) << total_mass_all() << "\n";
        }

        // ✅ CORRECTED: Final results show all masses
        void print_final_results() const {
            std::cout << "\n======================\n";
            std::cout << "FINAL MASSES PER FOLLICLE:\n";
            auto ms = masses();
            for (int f = 0; f < cfg_.Nf; ++f) {
                std::cout << "Follicle " << f << ": " << std::setw(10) << ms[f] << "\n";
            }
            std::cout << "TOTAL MASS: " << total_mass_all() << "\n";
            std::cout << "CSV files ready for Python visualization\n";
        }
};

int main() {
    FvmSolver::Config cfg;
    cfg.Nm = 2000;
    cfg.Nf = 3;
    cfg.T = 2.0;
    cfg.CFL = 0.3;
    cfg.output_freq = 20;

    FvmSolver solver(cfg);
    solver.solve_and_output();

    return 0;
}

