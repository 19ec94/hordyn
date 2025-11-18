#include "domain.hpp"
#include <stdexcept>
#include <cmath>

Domain::Domain(uint32_t num_cycles_, double cycle_duration_, double threshold_y_,
				uint32_t cells_per_half_cycle_, uint32_t num_cells_y_)
		: num_cycles(num_cycles_),
		cycle_duration(cycle_duration_),
		threshold_y(threshold_y_),
		cells_per_half_cycle(cells_per_half_cycle_),
		num_cells_y(num_cells_y_),
		num_cells_x(2 * cells_per_half_cycle_ * num_cycles_),
		length_x(cycle_duration_ * num_cycles_),
		length_y(1.0),
		dx(length_x / num_cells_x),
		dy(length_y / num_cells_y),
		threshold_y_index_y(std::round(threshold_y/dy)),
		coords_x(num_cells_x + 1),
		coords_y(num_cells_y + 1),
		centers_x(num_cells_x),
		centers_y(num_cells_y)
{
		if (num_cycles <= 0)
				throw std::invalid_argument("num_cycles must be > 0");
		if (cycle_duration <= 0)
				throw std::invalid_argument("cycle_duration must be > 0");
		if (threshold_y <= 0 || threshold_y > 1.0)
				throw std::invalid_argument("threshold_y must be in (0,1]");
		if (cells_per_half_cycle <= 0)
				throw std::invalid_argument("cells_per_half_cycle must be > 0");
		if (num_cells_y <= 0)
				throw std::invalid_argument("num_cells_y must be > 0");
		double y_s = threshold_y / dy;
		if (std::abs(y_s - std::round(y_s)) > 1e-12) {
				throw std::runtime_error("Ys must align exactly with vertical grid lines");
		}
		if (num_cells_x * num_cells_y > 1e8) {
				throw std::runtime_error("Grid size too large for memory");
		}

		initialize_grid();
}


void Domain::initialize_grid() {
		for (uint32_t i = 0; i <= num_cells_x; ++i) {
				coords_x[i] = i * dx;
		}
		for (uint32_t j = 0; j <= num_cells_y; ++j) {
				coords_y[j] = j * dy;
		}
		for (uint32_t i = 0; i < num_cells_x; ++i) {
				centers_x[i] = 0.5 * (coords_x[i] + coords_x[i + 1]);
		}
		for (uint32_t j = 0; j < num_cells_y; ++j) {
				centers_y[j] = 0.5 * (coords_y[j] + coords_y[j + 1]);
		}
}


uint32_t Domain::cell_cycle_index(double x) const {
		int p = static_cast<int>(x / length_x * num_cycles);
		if (p < 0)
				p = 0;
		if (static_cast<uint32_t>(p) >= num_cycles)
				p = num_cycles - 1;
		return p;
}

uint32_t Domain::get_num_cols() const {
		return num_cells_x;
}

uint32_t Domain::get_num_rows() const {
		return num_cells_y;
}

std::vector<uint32_t> Domain::interface_x_index_continuity() const {
		std::vector<uint32_t> indices(num_cycles, 0.0);
		for (uint32_t c = 0; c < num_cycles; ++c) {
				indices[c] = c  * 2 * cells_per_half_cycle  + cells_per_half_cycle;
		}
		return indices;
}

std::vector<uint32_t> Domain::interface_x_index_doubling() const {
		std::vector<uint32_t> indices(num_cycles, 0.0);
		for (uint32_t c = 0; c < num_cycles; ++c) {
				indices[c] = (c + 1 ) * 2 * cells_per_half_cycle;
		}
		return indices;
}

std::vector<uint32_t> Domain::cell_x_index_dirichlet() const {
		std::vector<uint32_t> indices(num_cycles * cells_per_half_cycle, 0.0);
		for (uint32_t c = 0; c < num_cycles; ++c) {
				for (uint32_t idx = 0; idx < cells_per_half_cycle; idx++){
						indices[ c * cells_per_half_cycle + idx] = c * 2 * cells_per_half_cycle + cells_per_half_cycle + idx;
				}
		}
		return indices;
}

