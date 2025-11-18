#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include <cstdint>

class Domain {
		public:
				// User input parameters (const)
				const uint32_t num_cycles;
				const double cycle_duration;
				const double threshold_y;
				const uint32_t cells_per_half_cycle;
				const uint32_t num_cells_y;
				// End of User input

				// Computed parameters
				const uint32_t num_cells_x;
				const double length_x;
				const double length_y;
				const double dx;
				const double dy;
				const uint32_t threshold_y_index_y;

				// Grid data
				std::vector<double> coords_x;
				std::vector<double> coords_y;
				std::vector<double> centers_x;
				std::vector<double> centers_y;

				// Constructor declaration
				Domain(uint32_t num_cycles_, double cycle_duration_, double threshold_y_,
								uint32_t cells_per_half_cycle_, uint32_t num_cells_y_);

				// Member function declarations
				void initialize_grid();
				uint32_t cell_cycle_index(double x) const;
				uint32_t get_num_cols() const;
				uint32_t get_num_rows() const;
				std::vector<uint32_t> interface_x_index_continuity() const;
				std::vector<uint32_t> interface_x_index_doubling() const;
				std::vector<uint32_t> cell_x_index_dirichlet() const;
};

#endif // DOMAIN_HPP
