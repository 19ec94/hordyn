#include "computations.hpp"

void compute_phi(const Domain& domain, std::vector<double>& phi,
        const FollicleParams& follicle_params) {
    const double cx = follicle_params.cx;
    const double cy = follicle_params.cy;
    const double sigma  = follicle_params.sigma;
    // Identify only Omega1 in the first the cycle (cycle starts from 0)
    const double factor = 1/(2 * std::numbers::pi * sigma * sigma);
    const uint32_t omega1_end = domain.cells_per_half_cycle;
    const uint32_t ys_idx = domain.ys_idx;
    const uint32_t num_cols = domain.Nx;

    for (uint32_t row = 0; row < ys_idx ; ++row) {
        double dy = domain.centers_y[row] - cy;
        for (uint32_t col = 0; col < omega1_end; ++col) {
            double dx = domain.centers_x[col] - cx;	
            size_t index = row * num_cols + col;
            phi[index] = factor * std::exp(-(dx*dx + dy*dy)/(2* sigma * sigma));
        }
    }
}

double compute_mass(const Domain& domain, const std::vector<double>& phi) {
    double result = 0.0;
    const uint32_t num_cols = domain.Nx;
    const uint32_t num_rows = domain.Ny;
    const double dx = domain.dx;
    const double dy = domain.dy;
    for (uint32_t i = 0; i < num_rows * num_cols; ++i) {
        result += phi[i] * dx * dy;
    }
    return result;
}


double compute_mf(const Domain& domain, const std::vector<double>& phi) {
    // TODO: Simplify
    double result = 0.0;
    const uint32_t num_cols = domain.Nx;
    const uint32_t num_rows = domain.Ny;
    const double dx = domain.dx;
    const double dy = domain.dy;
    for (uint32_t row = 0; row < num_rows ; ++row) {
        for (uint32_t col = 0; col < num_cols; ++col) {
            result += domain.centers_y[row] * phi[row * num_cols + col] * dx * dy;
        }
    }
    return result;
}


double compute_U(const GlobalBioParams& global_bio_params, const double M){
    double result = 0.0;
    const double numerator = 1.0 - global_bio_params.U_min; 
    const double denominator = 1.0 + std::exp(global_bio_params.c * (
                M - global_bio_params.M_ref)); 
    result = global_bio_params.U_min + (numerator / denominator);
    return result;
}


double compute_uf(const FollicleParams& follicle_params, const double mf, const double U){
    double result = follicle_params.b1 + 
        (std::exp(follicle_params.b2 * mf) / follicle_params.b3);
    if (result > 1.0) result = 1.0;
    return (result * U);
}

void compute_g_h(const Domain& domain,
        std::vector<double>& g, std::vector<double>& h,
        const FollicleParams& follicle_params,
        const double uf) {
    const uint32_t num_rows = domain.Ny;
    const uint32_t num_cols = domain.Nx;
    const uint32_t cells_per_half_cycle = domain.cells_per_half_cycle;
    const double ys_idx = domain.ys_idx;

    // Parameters for age veclocity
    const double gamma1 = follicle_params.gamma1;
    const double gamma2 = follicle_params.gamma2;

    // Parameters for maturation velocity
    const double c1 = follicle_params.c1;
    const double c2 = follicle_params.c2;
    const double tau_h = follicle_params.tau_h;
    const double u_bar = follicle_params.u_bar;

    for (uint32_t row = 0; row < num_rows; ++row) {
        double y = domain.centers_y[row];
        for (uint32_t col = 0; col < num_cols; ++col) {
            size_t index = row * num_cols + col;
            if (row < ys_idx) {
                if ( col < cells_per_half_cycle){
                    // Omega1
                    g[index] = gamma1 * uf + gamma2;
                    h[index] = tau_h * (-y * y + (c1 * y + c2) * (1.0 - std::exp(-uf / u_bar)));
                } else {
                    // Omega2
                    g[index] = 1.0;
                    h[index] = 0.0;
                } 
            } else {
                // Omega3
                g[index] = 1.0;
                h[index] = tau_h * (-y * y + (c1 * y + c2) * (1.0 - std::exp(-uf / u_bar)));
            }
        }
    }
}

void compute_source(const Domain& domain,
        std::vector<double>& source,
        const GlobalBioParams& global_bio_params, const double U) {
    const uint32_t num_rows = domain.Ny;
    const uint32_t num_cols = domain.Nx;
    const uint32_t cells_per_half_cycle = domain.cells_per_half_cycle;
    const double ys_idx = domain.ys_idx;

    const double Lambda_bar = global_bio_params.Lambda_bar;
    const double gamma_bar = global_bio_params.gamma_bar;
    const double ys = domain.ys;

    for (uint32_t row = 0; row < num_rows; ++row) {
        double y = domain.centers_y[row];
        for (uint32_t col = 0; col < num_cols; ++col) {
            size_t index = row * num_cols + col;
            if (row < ys_idx) {
                if ( col < cells_per_half_cycle){
                    // Omega1
                    source[index] = Lambda_bar * std::exp(-(pow(y - ys, 2) / pow(gamma_bar,2))) * (1 - U);
                } else {
                    // Omega2
                    source[index] = 0.0;
                } 
            } else {
                // Omega3
                source[index] = Lambda_bar * std::exp(-(pow(y - ys, 2) / pow(gamma_bar,2))) * (1 - U);
            }
        }
    }
}
