#include "computations.hpp"

#include <cmath>
#include <numbers>


/*
void compute_phi(const Domain& domain,
			   	std::vector<double>& phi,
			   	double cx, double cy, double sigma) {
		// Identify only Omega1 in the first the cycle (cycle starts from 0)
		const double factor = 1/(2 * std::numbers::pi * sigma * sigma);
		const uint32_t cells_per_half_cycle = domain.cells_per_half_cycle;
		const uint32_t threshold_y_index_y = domain.threshold_y_index_y;
		const uint32_t num_cols = domain.get_num_cols();

		for (uint32_t row = 0; row < threshold_y_index_y; ++row) {
				//double dy = domain.centers_y[row] - cy;
				double value = 0.0;
				double y = domain.centers_y[row];
						if (y >= 0.05 && y < 0.1) {
								value = 8;
						} else if (y >= 0.1 && y < 0.15) {
								value = 7;
						} else if (y >= 0.15 && y < 0.2) {
								value = 5;
						} else {
								value = 0;
						}

				for (uint32_t col = 0; col < cells_per_half_cycle * 2 ; ++col) {
						//double dx = domain.centers_x[col] - cx;	
						size_t index = row * num_cols + col;
						phi[index] = value; 
				}
		}
}

*/

void compute_phi(const Domain& domain, std::vector<double>& phi, double cx, 
				double cy, double sigma) {
		// Identify only Omega1 in the first the cycle (cycle starts from 0)
		const double factor = 1/(2 * std::numbers::pi * sigma * sigma);
		const uint32_t cells_per_half_cycle = domain.cells_per_half_cycle;
		const uint32_t threshold_y_index_y = domain.threshold_y_index_y;
		const uint32_t num_cols = domain.get_num_cols();

		for (uint32_t row = 0; row < threshold_y_index_y; ++row) {
				double dy = domain.centers_y[row] - cy;
				for (uint32_t col = 0; col < cells_per_half_cycle ; ++col) {
						double dx = domain.centers_x[col] - cx;	
						size_t index = row * num_cols + col;
						phi[index] = factor * std::exp(-(dx*dx + dy*dy)/(2* sigma * sigma));
				}
		}
}

double compute_mass(const Domain& domain,
			   	const std::vector<double>& phi) {
		double result = 0.0;
		const uint32_t num_cols = domain.get_num_cols();
		const uint32_t num_rows = domain.get_num_rows();
		const double dx = domain.dx;
		const double dy = domain.dy;
		for (uint32_t i = 0; i < num_rows * num_cols; ++i) {
				result += phi[i] * dx * dy;
		}
		return result;
}


double compute_local_maturity(const Domain& domain,
			   	const std::vector<double>& phi) {
		double result = 0.0;
		const uint32_t num_cols = domain.get_num_cols();
		const uint32_t num_rows = domain.get_num_rows();
		const double dx = domain.dx;
		const double dy = domain.dy;
		for (uint32_t row = 0; row < num_rows ; ++row) {
				for (uint32_t col = 0; col < num_cols; ++col) {
						result += domain.centers_y[row] * phi[row * num_cols + col] * dx * dy;
				}
		}
		return result;
}


double compute_global_fsh(const GlobalBioParams& global_bio_params,
			   	const double global_maturity){
		double result = 0.0;
		const double numerator = 1.0 - global_bio_params.U_min; 
		const double denominator = 1.0 + std::exp(global_bio_params.c * (
								global_maturity - global_bio_params.M_ref)); 
		result = global_bio_params.U_min + (numerator / denominator);
		return result;
}


double compute_local_fsh(const FollicleParams& follicle_params,
			   	const double local_maturity,
				const double global_fsh){
		double result = follicle_params.b1 + 
			   	(std::exp(follicle_params.b2 * local_maturity) / follicle_params.b3);
		if (result > 1.0) result = 1.0;
		return (result * global_fsh);
}

void compute_g_h_source(const Domain& domain,
			   	std::vector<double>& g,
			   	std::vector<double>& h,
				std::vector<double>& source,
			   	const GlobalBioParams& global_bio_params, 
				const FollicleParams& follicle_params,
				const double local_fsh,
			   	const double global_fsh) {
		const uint32_t num_rows = domain.get_num_rows();
		const uint32_t num_cols = domain.get_num_cols();
		const double cd = domain.cycle_duration;
		for (uint32_t row = 0; row < num_rows; ++row) {
				double y = domain.centers_y[row];
				for (uint32_t col = 0; col < num_cols; ++col) {
						double x = domain.centers_x[col];
						uint32_t cycle = domain.cell_cycle_index(x);
						size_t index = row * num_cols + col;
						if (y < domain.threshold_y) {
								if ( (x >=  cycle * cd) && (x < (cycle + 0.5) * cd)){
										// Omega1
										//g[index] = 1.0;
										//h[index] = 0.0;
										//source[index] = 0.0;
										g[index] = follicle_params.gamma1 * local_fsh + follicle_params.gamma2;
										h[index] = follicle_params.tau_h * (-y * y + (follicle_params.c1 * y + follicle_params.c2) * (1.0 - std::exp(-local_fsh / follicle_params.u_bar)));
										source[index] = global_bio_params.Lambda_bar * std::exp(-(pow(y - global_bio_params.y_s, 2) / pow(global_bio_params.gamma_bar,2))) * (1 - global_fsh);
								} else if ( (x >= (cycle + 0.5) * cd )  && (x < (cycle+1) * cd )) {
										// Omega2
										g[index] = 1.0;
										h[index] = 0.0;
										source[index] = 0.0;
								} else {
										// Edge case
								}
						} else {
								// Omega3
								//g[index] = 0.0;
								//h[index] = 0.0;
								//source[index] = 0.0;
								g[index] = 1.0;
								h[index] = follicle_params.tau_h * (-y * y + (follicle_params.c1 * y + follicle_params.c2) * (1.0 - std::exp(-local_fsh / follicle_params.u_bar)));
								source[index] = global_bio_params.Lambda_bar * std::exp(-(pow(y - global_bio_params.y_s, 2) / pow(global_bio_params.gamma_bar,2))) * (1 - global_fsh);
						}
				}
		}
}
