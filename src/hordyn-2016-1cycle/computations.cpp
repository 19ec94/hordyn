#include "computations.hpp"
#include <iostream>

void compute_phi(const Domain& domain, std::vector<double>& phi,
        const FollicleParams& follicle_params) {

    const double M0 = 1.0;
    //const double k = log(2)/domain.Dc;
    double g = 0.5;
    const double CG1 = (2 * log(2)) / (domain.Dc * (2 + pow(2, (1-0.5)/domain.Dc) *(g-1.0) -g));
    const double CSM = g * CG1;

    const uint32_t num_rows = domain.Ny;
    const uint32_t num_cols = domain.Nx;
    const uint32_t cells_per_half_cycle = domain.cells_per_half_cycle;
    const double ys_idx = domain.ys_idx;

    const double mu1 = follicle_params.mu1;
    const double mu2 = follicle_params.mu2;

    const int s1 = static_cast<int>(std::round(mu1/domain.dy));
    const int s2 = static_cast<int>(std::round(mu2/domain.dy));

    for (uint32_t row = 0; row < num_rows; ++row) {
        for (uint32_t col = 0; col < num_cols; ++col) {
            double x = domain.centers_x[col];
            size_t index = row * num_cols + col;
            //if ( row < ys_idx) {
            if (row >= s1 && row < s2) {
                if ( col < cells_per_half_cycle){
                    // Omega1
                    phi[index] = (1.0/(mu2-mu1)) * (M0 * CG1 * pow(2, -x/domain.Dc));
                    //phi[index] = 1.0;
                } else {
                    // Omega2
                    phi[index] = (1.0/(mu2-mu1)) * (M0 * CSM * pow(2, -x/domain.Dc));
                    //phi[index] = 0.5;
                } 
            }
        }
    }
}

void compute_phi_old(const Domain& domain, std::vector<double>& phi,
        const FollicleParams& follicle_params) {

    const double M0 = 1.0;
    //const double k = log(2)/domain.Dc;
    double g = 0.5;
    const double CG1 = (2 * log(2)) / (domain.Dc * (2 + pow(2, (1-0.5)/domain.Dc) *(g-1.0) -g));
    const double CSM = g * CG1;

    const uint32_t num_rows = domain.Ny;
    const uint32_t num_cols = domain.Nx;
    const uint32_t cells_per_half_cycle = domain.cells_per_half_cycle;
    const double ys_idx = domain.ys_idx;

    for (uint32_t row = 0; row < num_rows; ++row) {
        for (uint32_t col = 0; col < num_cols; ++col) {
            double x = domain.centers_x[col];
            size_t index = row * num_cols + col;
            if (row < ys_idx) {
                if ( col < cells_per_half_cycle){
                    // Omega1
                    phi[index] = M0 * CG1 * pow(2, -x/domain.Dc);
                    //phi[index] = 1.0;
                } else {
                    // Omega2
                    phi[index] = M0 * CSM * pow(2, -x/domain.Dc);
                    //phi[index] = 0.5;
                } 
            }
        }
    }
}

/*
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
*/

double compute_mass(const Domain& domain, const std::vector<double>& phi) {
    double result = 0.0;
    const uint32_t num_cols = domain.Nx;
    const uint32_t num_rows = domain.Ny;
    const double dx = domain.dx;
    const double dy = domain.dy;
    for (uint32_t i = 0; i < num_rows * num_cols; ++i) {
        result += phi[i];
    }
    result *= dx * dy;
    return result;
}

double compute_mass_mitosis(const Domain& domain, const std::vector<double>& phi){
    double result = 0.0;
    const uint32_t num_cols = domain.Nx;
    const uint32_t num_rows = domain.Ny;
    const double dx = domain.dx;
    const double dy = domain.dy;
    const uint32_t ys_idx = domain.ys_idx;
    for (uint32_t i = 0; i < ys_idx; i++) {
        // TODO adjust the contribution to match 30 minutes
        for (uint32_t j = num_cols-1;  j > num_cols-2; j--) {
            uint32_t index = i * num_cols + j;
            result += phi[index] * dx * dy;
        }
    }
    return result;
}

double compute_mass_omega12(const Domain& domain, const std::vector<double>& phi){
    double result = 0.0;
    const uint32_t num_cols = domain.Nx;
    const uint32_t num_rows = domain.Ny;
    const double dx = domain.dx;
    const double dy = domain.dy;
    const uint32_t ys_idx = domain.ys_idx;
    for (uint32_t i = 0; i < ys_idx; i++) {
        for (uint32_t j = 0; j < num_cols; j++) {
            uint32_t index = i * num_cols + j;
            result += phi[index] * dx * dy;
        }
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

double compute_U_ol(const GlobalBioParams& global_bio_params, double ts, double t){
    double result = 0.0;
    if (t <= ts) result = global_bio_params.U_max;
    else if (t > ts) result = global_bio_params.U_min;
    return result;
}


double compute_uf(const FollicleParams& follicle_params, const double mf, const double U){
    double result = follicle_params.b1 + 
        (std::exp(follicle_params.b2 * mf) / follicle_params.b3);
    if (result > 1.0) result = 1.0;
    return (result * U);
}

double compute_uf_ol(const GlobalBioParams& global_bio_params, double ts, double t){
    double result = 0.0;
    if (t <= ts) result = global_bio_params.U_max;
    else if (t > ts) result = global_bio_params.U_min;
    return result;
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
                    //g[index] = gamma1 * uf + gamma2;
                    h[index] = tau_h * (-y * y + (c1 * y + c2) * (1.0 - std::exp(-uf / u_bar)));
                    g[index] = 0.5;
                    //h[index] = 0.0;
                } else {
                    // Omega2
                    g[index] = 1.0;
                    h[index] = 0.0;
                } 
            } else {
                // Omega3
                //g[index] = 1.0;
                h[index] = tau_h * (-y * y + (c1 * y + c2) * (1.0 - std::exp(-uf / u_bar)));
                g[index] = 0.0;
                //h[index] = 0.0;
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

    // TODO: Make it general
    const int ny = static_cast<int>(std::round((0.55 - domain.ys)/domain.dy));
    //std::cout << "ny is " << ny << "\n";

    const double U_max = global_bio_params.U_max;

    for (uint32_t row = 0; row < num_rows; ++row) {
        double y = domain.centers_y[row];
        for (uint32_t col = 0; col < num_cols; ++col) {
            size_t index = row * num_cols + col;
            if ((row >= ys_idx - ny) && (row < ys_idx)) {
                if ( col < cells_per_half_cycle){
                    // Omega1
                    //source[index] = Lambda_bar * std::exp(-(pow(y - ys, 2) / pow(gamma_bar,2))) * ((U_max - U)/U_max);
                    source[index] = 0.0;
                } else {
                    // Omega2
                    source[index] = 0.0;
                } 
            } else if ( (row >= ys_idx) && (row <  ys_idx + ny)){
                // Omega3
                //source[index] = Lambda_bar * std::exp(-(pow(y - ys, 2) / pow(gamma_bar,2))) * ((U_max - U)/U_max);
                source[index] = 0.0;
            } else {
                source[index] = 0.0;
            }
        }
    }
}
