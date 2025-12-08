#pragma once

#include <vector>
#include <cstdint>

class Domain {
		public:
				// User input parameters (const)
				const double Dc;
				const double ys;
				const uint32_t cells_per_half_cycle;
				const uint32_t Ny;
				// End of user input
				
				// Compute parameters
				const uint32_t Nx;
				const double Lx, Ly;
				const double dx, dy;
				const uint32_t ys_idx;

				// Grid data
				std::vector<double> coords_x, coords_y;
				std::vector<double> centers_x, centers_y;

				// Constructor declaration
				Domain(double Dc_, double ys_, uint32_t cells_per_half_cycle_, uint32_t Ny_);

				void initialize_grid();
};
