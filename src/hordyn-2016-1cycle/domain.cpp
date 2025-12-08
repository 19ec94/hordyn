#include "domain.hpp"
#include <stdexcept>
#include <cmath>


Domain::Domain(double Dc_, double ys_, uint32_t cells_per_half_cycle_, uint32_t Ny_)
    : Dc(Dc_), ys(ys_), cells_per_half_cycle(cells_per_half_cycle_), Ny(Ny_),
    Nx(2 * cells_per_half_cycle),
    Lx(Dc_), Ly(1.0), 
    dx(Lx / Nx), dy(Ly / Ny),
    ys_idx(std::round(ys/dy)),
    coords_x(Nx + 1), coords_y(Ny +1), centers_x(Nx), centers_y(Ny) {
        if (Dc <= 0)
            throw std::invalid_argument("Dc must be >0.");
        if (ys <= 0 || ys > 1)
            throw std::invalid_argument("ys must be in ]0, 1[.");
        if (Nx <=0)
            throw std::invalid_argument("Nx must be > 0");
        if (Ny <= 0)
            throw std::invalid_argument("Ny must be >0");

        double ys_calculated = ys / dy;
        if (std::abs(ys_calculated - std::round(ys_calculated)) > 1e-12) 
            throw std::runtime_error("ys must along with the grid line.");

        if (Nx * Ny > 1e8) 
            throw std::runtime_error("Grid size too large for memory");

        initialize_grid();
    }

void Domain::initialize_grid(){
    for (uint32_t i = 0; i < Nx+1; ++i) coords_x[i] = i * dx;
    for (uint32_t i = 0; i < Nx; ++i) centers_x[i] = (i * dx) + (dx/2);
    for (uint32_t j = 0; j < Ny+1; ++j) coords_y[j] = j * dy;
    for (uint32_t j = 0; j < Ny; ++j) centers_y[j] = (j * dy) + (dy/2);
}

