#include "compute_rhs.hpp"
#include <algorithm>
#include <cmath>


std::vector<double> compute_rhs(const Domain& domain,
				const Time& time,
				const std::vector<double>& phi,
				const std::vector<double>& g,
				const std::vector<double>& h,
				const std::vector<double>& source
				) {
		double dx = domain.dx;
		double dy = domain.dy;
		double dt = time.dt;
		int num_rows = domain.get_num_rows();  // Number of rows
		int num_cols = domain.get_num_cols();  // Number of columns
		int interfaceYIndexDirichlet = domain.threshold_y_index_y;
		std::vector<uint32_t> interface_x_index_doubling = domain.interface_x_index_doubling();
		std::vector<uint32_t> cell_x_index_dirichlet = domain.cell_x_index_dirichlet();
		// Index helpers
		auto idx = [num_cols](int i, int j) { return i * num_cols + j; };
		// Vertical interface fluxes size num_rows*(num_cols+1)
		auto idx_G = [num_cols](int i, int j) { return i * (num_cols + 1) + j; };
		// Horizontal interface fluxes size (num_rows+1)*num_cols flattened
		//auto idx_H = [num_cols](int i, int j) { return i * num_cols + j; };      

		// Flux arrays for vertical and horizontal interfaces
		std::vector<double> G_left(num_rows * (num_cols + 1), 0.0);
		std::vector<double> G_right(num_rows * (num_cols + 1), 0.0);
		std::vector<double> H_top((num_rows + 1) * num_cols, 0.0);
		std::vector<double> H_bottom((num_rows + 1) * num_cols, 0.0);

		// Function to compute slope ratio for limiters
		auto compute_rk = [&](double vel0, double vel1, double vel2, double vel3,
						double phi0, double phi1, double phi2, double phi3) -> double {
				if (vel0 >= 0 && vel1 >= 0 && vel2 >= 0) {
						double num = (vel1 * phi1) - (vel0 * phi0);
						double den = (vel2 * phi2) - (vel1 * phi1);
						if (std::fabs(den) < 1e-14) return 0.0;
						return num / den;
				} else if (vel1 <= 0 && vel2 <= 0 && vel3 <= 0) {
						double num = (vel3 * phi3) - (vel2 * phi2);
						double den = (vel2 * phi2) - (vel1 * phi1);
						if (std::fabs(den) < 1e-14) return 0.0;
						return num / den;
				} else {
						return 0.0;
				}
		};

		// Korean limiter function
		auto compute_korean_limiter = [](double slope) -> double {
				double intermediate = std::min({2.0 * slope, (2.0 + slope) / 3.0, 2.0});
				return std::max(0.0, intermediate);
		};

		// Compute vertical fluxes across interfaces
		for (int row = 0; row < num_rows; ++row) {
				std::vector<bool> isFluxDoubling(num_cols + 1, false);
				for (int col = 0; col <= num_cols; ++col) {
						if ((row < interfaceYIndexDirichlet) &&
										(std::find(interface_x_index_doubling.begin(), interface_x_index_doubling.end(), col) != interface_x_index_doubling.end())) {
								isFluxDoubling[col] = true;
						}
				}

				for (int col = 1; col < num_cols; ++col) {
						double g_low = 0.0, g_high = 0.0, rk = 0.0, limiter = 0.0;

						if (isFluxDoubling[col]) {
								g_low = std::max(g[idx(row, col - 1)], 0.0) * phi[idx(row, col - 1)]
										+ (std::min(g[idx(row, col)], 0.0) * phi[idx(row, col)]) / 2.0;
								g_high = (g[idx(row, col - 1)] * phi[idx(row, col - 1)] + (g[idx(row, col)] * phi[idx(row, col)] / 2.0)) / 2.0;
						} else {
								g_low = std::max(g[idx(row, col - 1)], 0.0) * phi[idx(row, col - 1)]
										+ std::min(g[idx(row, col)], 0.0) * phi[idx(row, col)];
								g_high = (g[idx(row, col - 1)] * phi[idx(row, col - 1)] + g[idx(row, col)] * phi[idx(row, col)]) / 2.0;
						}

						if (col >= 2 && col <= num_cols - 2) {
								int K = -1;
								for (int k = col - 1; k <= col + 1; ++k) {
										if (k >= 1 && k < num_cols && isFluxDoubling[k]) {
												K = k;
												break;
										}
								}
								if (K == -1) {
										rk = compute_rk(
														g[idx(row, col - 2)], g[idx(row, col - 1)], g[idx(row, col)], g[idx(row, col + 1)],
														phi[idx(row, col - 2)], phi[idx(row, col - 1)], phi[idx(row, col)], phi[idx(row, col + 1)]
													   );
								} else {
										if (col == K - 1) {
												rk = compute_rk(
																g[idx(row, K - 3)], g[idx(row, K - 2)], g[idx(row, K - 1)], 0.5 * g[idx(row, K)],
																phi[idx(row, K - 3)], phi[idx(row, K - 2)], phi[idx(row, K - 1)], phi[idx(row, K)]
															   );
										} else if (col == K) {
												rk = compute_rk(
																2.0 * g[idx(row, K - 2)], 2.0 * g[idx(row, K - 1)], g[idx(row, K)], g[idx(row, K + 1)],
																phi[idx(row, K - 2)], phi[idx(row, K - 1)], phi[idx(row, K)], phi[idx(row, K + 1)]
															   );
										} else if (col == K + 1) {
												rk = compute_rk(
																2.0 * g[idx(row, K - 1)], g[idx(row, K)], g[idx(row, K + 1)], g[idx(row, K + 2)],
																phi[idx(row, K - 1)], phi[idx(row, K)], phi[idx(row, K + 1)], phi[idx(row, K + 2)]
															   );
										} else {
												// TODO: Throw runtime error
												rk = compute_rk(
																g[idx(row, col - 2)], g[idx(row, col - 1)], g[idx(row, col)], g[idx(row, col + 1)],
																phi[idx(row, col - 2)], phi[idx(row, col - 1)], phi[idx(row, col)], phi[idx(row, col + 1)]
															   );
										}
								}
						} if (col == num_cols-1) {
										rk = compute_rk(
														g[idx(row, num_cols - 3)], g[idx(row, num_cols - 2)], g[idx(row, num_cols - 1)], 0.5 * g[idx(row, 0)],
														phi[idx(row, num_cols - 3)], phi[idx(row, num_cols - 2)], phi[idx(row, num_cols - 1)], phi[idx(row, 0)]
													   );
						} else if (col == 1) {
										rk = compute_rk(
														2.0 * g[idx(row, num_cols - 1)], g[idx(row, 0)], g[idx(row, 1)], g[idx(row, 2)],
														phi[idx(row, num_cols - 1)], phi[idx(row, 0)], phi[idx(row, 1)], phi[idx(row, 2)]
													   );
						}

						limiter = compute_korean_limiter(rk);

						double flux_left = g_low + limiter * (g_high - g_low);
						G_left[idx_G(row, col)] = flux_left;
						G_right[idx_G(row, col)] = isFluxDoubling[col] ? 2.0 * flux_left : flux_left;
				}
		}
		// Periodic BC in the age axis 
		for (int row = 0; row < num_rows; ++row) {
				double g_low = std::max(g[idx(row, num_cols - 1)], 0.0) * phi[idx(row, num_cols - 1)]
						+ (std::min(g[idx(row, 0)], 0.0) * phi[idx(row, 0)]) / 2.0;
				double g_high = (g[idx(row, num_cols - 1)] * phi[idx(row, num_cols - 1)] + (g[idx(row, 0)] * phi[idx(row, 0)] / 2.0)) / 2.0;
				double rk = compute_rk(
								2.0 * g[idx(row, num_cols - 2)], 2.0 * g[idx(row, num_cols - 1)], g[idx(row, 0)], g[idx(row, 1)],
								phi[idx(row, num_cols - 2)], phi[idx(row, num_cols - 1)], phi[idx(row, 0)], phi[idx(row, 1)]
							   );
				double limiter = compute_korean_limiter(rk);
				double flux_left = g_low + limiter * (g_high - g_low);
				G_left[idx_G(row, num_cols)] = flux_left;
				G_right[idx_G(row, num_cols)] = 2.0 * flux_left;
				// data at first interface  ==  data at last interface
				G_left[idx_G(row, 0)] = G_left[idx_G(row, num_cols)];
				G_right[idx_G(row, 0)] = G_right[idx_G(row, num_cols)];
		}

		// Compute horizontal fluxes across interfaces
		for (int col = 0; col < num_cols; ++col) {
				for (int row = 1; row < num_rows; ++row) {
						double h_low = 0.0, h_high = 0.0, rk = 0.0, limiter = 0.0;

						if ((row == interfaceYIndexDirichlet) &&
										(std::find(cell_x_index_dirichlet.begin(), cell_x_index_dirichlet.end(), col) != cell_x_index_dirichlet.end())) {
								h_low = 0.0;
								h_high = 0.0;
						} else {
								h_low = std::max(h[idx(row - 1, col)], 0.0) * phi[idx(row - 1, col)]
										+ std::min(h[idx(row, col)], 0.0) * phi[idx(row, col)];
								h_high = (h[idx(row - 1, col)] * phi[idx(row - 1, col)] + h[idx(row, col)] * phi[idx(row, col)]) / 2.0;
						}

						if (row >= 2 && row <= num_rows - 2) {
								rk = compute_rk(
												h[idx(row - 2, col)], h[idx(row - 1, col)], h[idx(row, col)], h[idx(row + 1, col)],
												phi[idx(row - 2, col)], phi[idx(row - 1, col)], phi[idx(row, col)], phi[idx(row + 1, col)]
											   );
						}

						limiter = compute_korean_limiter(rk);

						double flux_top = h_low + limiter * (h_high - h_low);

						H_top[row * num_cols + col] = flux_top;
						H_bottom[row * num_cols + col] = flux_top;
				}
		}

		// Compute RHS values for each cell using flux differences and source
		std::vector<double> rhs(num_rows * num_cols, 0.0);
		for (int row = 0; row < num_rows; ++row) {
				for (int col = 0; col < num_cols; ++col) {
						rhs[idx(row, col)] =
								(dt / dx) * (G_left[idx_G(row, col + 1)] - G_right[idx_G(row, col)])
								+ (dt / dy) * (H_bottom[(row + 1) * num_cols + col] - H_top[row * num_cols + col])
								+ dt * source[idx(row, col)] * phi[idx(row, col)];
				}
		}

		return rhs;
}
