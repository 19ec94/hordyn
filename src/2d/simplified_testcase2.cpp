#include <iostream>
#include <vector>
#include <ranges>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <cstdio>
#include <numbers>

using GridField2D = std::vector<std::vector<double>>;

template<typename T>
void printMatrix(const std::string& name, const std::vector<std::vector<T>>& matrix) {
  std::cout << "******** " << name << " *******"  <<"\n";
  int numRows = matrix.size();
  int numCols = matrix[0].size();

  for (int c = numCols-1; c >= 0; --c){
    for (int r = 0; r < numRows; ++r) {
      std::cout << matrix[r][c] << " ";
    }
    std::cout << "\n";
  }
}

template<typename T>
void printVector(const std::string& name, const std::vector<T>& vec) {
  std::cout << "******** " << name << " *******"  <<"\n";
  for (const auto& elem: vec) {
    std::cout << elem << " ";
  }
  std::cout << "\n";
}


double evaluateGaussianQuadrature2D(double x0, double x1, double y0, double y1,
    double cx, double cy, double sigma)
{
  const double gaussPts[2] = { -0.5773502692, 0.5773502692 };
  const double gaussWts[2] = { 1.0, 1.0 };

  double dx = x1 - x0;
  double dy = y1 - y0;

  double integral = 0.0;
  // 2x2 Gaussian quadrature
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      double xq = 0.5 * (dx * gaussPts[i] + (x1 + x0));
      double yq = 0.5 * (dy * gaussPts[j] + (y1 + y0));

      double dx_center = xq - cx;
      double dy_center = yq - cy;

      double val = std::exp(-(dx_center*dx_center + dy_center*dy_center) / (2.0 * sigma * sigma));

      integral += gaussWts[i] * gaussWts[j] * val;
    }
  }

  integral *= 0.25 * dx * dy;  // Jacobian of transformation
  return integral / (dx * dy); // Average value over the cell
}


double findAbsoluteMaximumInMatrix(const GridField2D& matrix) {
  // Initialize the maximum value to  negative infinity
  double maxVal = -std::numeric_limits<double>::infinity(); 
  for (const auto& row: matrix) {
    if (!row.empty()) {
      // Find element with maximum absolute value in the row
      double rowMax = *std::max_element(row.begin(), row.end(), 
          [](double a, double b){ return std::abs(a) < std::abs(b); }
          );
      maxVal = (std::abs(rowMax) > maxVal) ? std::abs(rowMax) : maxVal;
    }
  }
  return maxVal;
}


GridField2D compute_rhs(int Nx, int Ny, 
    const GridField2D& phi,
    const GridField2D& g, 
    const GridField2D& h, 
    const GridField2D& Lambda,
    const std::vector<int>& interfaceXIndexContinuity,
    const std::vector<int>& interfaceXIndexDoubling,
    const std::vector<int>& cellXIndexDirichlet,
    const int interfaceYIndexDirichlet,
    const double dx,
    const double dy,
    const double dt
    ) {
  GridField2D G_left(Nx+1, std::vector<double>(Ny, 0.0));
  GridField2D G_right(Nx+1, std::vector<double>(Ny, 0.0));
  GridField2D H_top(Nx, std::vector<double>(Ny+1, 0.0));
  GridField2D H_bottom(Nx, std::vector<double>(Ny+1, 0.0));

  //auto periodicXIndex = [](int k, int Nx){return (k+Nx) % (Nx);};
  //auto periodicYIndex = [](int k, int Ny){return (k+Ny) % (Ny);};

  auto compute_rk = [](
      double vel0, double vel1, double vel2, double vel3, double phi0,
      double phi1, double phi2, double phi3){
    // TODO: Should vel3 be greater than or equal to zero? Refer the paper.
    if (vel0 >= 0 && vel1 >= 0 && vel2 >= 0) {
      double num = ( vel1 * phi1 ) - ( vel0 * phi0 );
      double den = ( vel2 * phi2 ) - ( vel1 * phi1 ) ;
      if (std::fabs(den) < 1e-14) { return 0.0; }
      return num / den;
      // TODO: Should vel0 be less than or equal to zero? Refer the paper.
    } else if (vel1 <=0 && vel2 <=0 && vel3 <= 0) {
      double num = ( vel3 * phi3 ) - ( vel2 * phi2 );
      double den = ( vel2 * phi2 ) - ( vel1 * phi1 );
      if (std::fabs(den) < 1e-14) { return 0.0; }
      return num / den;
    } else {
      return 0.0;
    }
  };
  auto compute_korean_limiter = [](double slope) { 
    double intermediate = std::min( {2.0 * slope, (2.0 + slope) / 3.0, 2.0 });
    return std::max(0.0, intermediate);
  };

  // 1. Compute flux across  Ny-1 vertical interfaces for each row
  for (int j = 0; j < Ny; ++j) {
    // Get which of the stencils (points) in the current row are FluxDoubling 
    // TODO: Check sizes
    std::vector<bool> isFluxDoubling(Nx + 1, false);
    //std::cout << j << "\n";
    for (int i = 0; i <= Nx; ++i) {
      if ( (j < interfaceYIndexDirichlet) && (std::find(interfaceXIndexDoubling.begin(), interfaceXIndexDoubling.end(), i) != interfaceXIndexDoubling.end()) ) {
        isFluxDoubling[i] = true;
      } else {
        isFluxDoubling[i] = false;
      }
    }

    // Exclude the flux at boundary interface 0 and Nx, treat them separately
    // below
    for (int  i = 1; i < Nx; ++i) {
      double g_low = 0.0;
      double g_high = 0.0;
      double rk = 0.0;
      double limiter = 0.0;
      // 1. Compute flux_low and flux_high
      if (isFluxDoubling[i]) {
        g_low  = std::max( g[i-1][j], 0.0) * phi[i-1][j]
          + (std::min(g[i][j], 0.0) * phi[i][j])/ 2.0;
        g_high  = (g[i-1][j] * phi[i-1][j] + ( (g[i][j] * phi[i][j])/2.0 ))/2.0;
      } else {
        g_low = std::max( g[i-1][j], 0.0) * phi[i-1][j] + std::min(g[i][j], 0.0) * phi[i][j];
        g_high = (g[i-1][j] * phi[i-1][j] + g[i][j] * phi[i][j])/2.0;
      }

      // 3. Calculate the limiter (excluding the boundary stenciles, and
      // a stencil close to the boundary. The limiter value there (e.g. at 0,
      // and 1 ) are zero, i.e r_0 = 0, r_1 = 0, r_{Nx-1} = 0, r_Nx = 0.0
      if (i >= 2 && i <= Nx-2) {
        // 3.1 Find flux doubling interface K near stencil i
        int K = -1; 
        for (int k = i-1; k <= i+1; ++k) {
          if (k >= 1 && k < Nx && isFluxDoubling[k] ) {
            K = k;
            break;
          }
        }
        if (K == -1) {
          // No FluxDoubling in stencil
          rk = compute_rk( 
              g[i-2][j], g[i-1][j], g[i][j], g[i+1][j], 
              phi[i-2][j], phi[i-1][j], phi[i][j], phi[i+1][j]
              );
        } else {
          if ( i == K-1 ) {
            rk = compute_rk( 
                g[K-3][j], g[K-2][j], g[K-1][j], 0.5 * g[K][j], 
                phi[K-3][j], phi[K-2][j], phi[K-1][j], phi[K][j]
                );
          } else if (i == K) {
            rk = compute_rk( 
                2.0 * g[K-2][j], 2.0 * g[K-1][j], g[K][j],  g[K+1][j], 
                phi[K-2][j], phi[K-1][j], phi[K][j], phi[K+1][j]
                );
          } else if ( i == K + 1) {
            rk = compute_rk( 
                2.0 * g[K-1][j], g[K][j], g[K+1][j],  g[K+2][j], 
                phi[K-1][j], phi[K][j], phi[K+1][j], phi[K+2][j]
                );
          } else {
            // TODO: Raise runtime error
            // fallback - no adjustment
            rk = compute_rk( 
                g[i-2][j], g[i-1][j], g[i][j], g[i+1][j], 
                phi[i-2][j], phi[i-1][j], phi[i][j], phi[i+1][j]
                );
          }
        }
      }
      // 3.2 Compute limiter
      limiter = compute_korean_limiter(rk);


      // 4. Assign left and right fluxes at interface i and j
      double flux_left = g_low + limiter * (g_high - g_low); 
      G_left[i][j] = flux_left; 
      G_right[i][j] = isFluxDoubling[i] ? 2.0 * flux_left : flux_left;
    }
  }

  //  Compute fluxes acorss Nx-1 row of  horizontal interfaces
  for (int i = 0; i < Nx; ++i) {
    double rk = 0.0;
    double limiter = 0.0;
    double h_low = 0.0;
    double h_high = 0.0;
    // compute flux_low and flux_high across the horizontal interface i,j
    for (int j = 1; j < Ny; ++j) {
      // if the interface is Dirichlet
      if ( (j == interfaceYIndexDirichlet) && (std::find(cellXIndexDirichlet.begin(),
              cellXIndexDirichlet.end(), i) !=
            cellXIndexDirichlet.end()) ) {
        h_low = 0.0;
        h_high = 0.0;
      } else {
        h_low = std::max(h[i][j-1], 0.0) * phi[i][j-1]
          + std::min(h[i][j], 0.0) * phi[i][j];
        h_high = (h[i][j-1] * phi[i][j-1] + h[i][j] * phi[i][j])/2.0;
      }

      // compute limiter excluding boundary interfaces and the one next to it
      if (j >= 2 && j <= Ny-2 ) {
        rk = compute_rk( h[i][j-2], h[i][j-1], h[i][j], h[i][j+1],
            phi[i][j-2], phi[i][j-1], phi[i][j], phi[i][j+1]);
      }
      limiter = compute_korean_limiter(rk);
      // 4. Assign top and bottom fluxes at interface i and j
      double flux_top = h_low + limiter * (h_high - h_low); 
      H_top[i][j] = flux_top;
      H_bottom[i][j] = flux_top;
    }
  }

  // 3. Compute rhs
  GridField2D rhs(Nx, std::vector<double>(Ny, 0.0));
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      rhs[i][j] = (dt/dx) * (G_left[i+1][j] - G_right[i][j])
        + (dt/dy) * (H_bottom[i][j+1] - H_top[i][j])
        + dt * Lambda[i][j] * phi[i][j];
    }
  }

  return rhs; 
}

void write_to_file(const std::string& filename,
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<double>>& phi
    ){
  // Output final phi at cell centers (no ghost)
  std::ofstream file(filename);
  file << "x,y,phi\n";
  for(int i=0; i < int(x.size()); ++i){
    for(int j=0; j < int(y.size()); ++j){
      file << std::setprecision(8) << x[i] << "," << y[j] << "," << phi[i][j] << "\n";
    }
  }
  file.close();
}

int main() {
  // Domain Parameters
  int Nc = 1; // Number of cell cycles
  double Dc = 1.0; 
  double Ys = 0.4;
  int cellsPerHalfCycle = 500;
  int Ny = 50;
  double Ly = 1.0;
  int Nx = 2 * Nc * cellsPerHalfCycle;
  double Lx = Dc * Nc;
  double dx = Lx / Nx;
  double dy = Ly / Ny;
  double j_s = Ys / dy;
  if (std::abs(j_s - std::round(j_s)) > 1e-12) {
    throw std::runtime_error("Ys must align exactly with vertical grid lines.");
  }
  if (Nx * Ny > 1e8) {
    throw std::runtime_error("Grid size too large for memory");
  }
  // Grid  
  std::vector<double> xCoords (Nx+1, 0.0);
  std::vector<double> yCoords (Ny+1, 0.0);
  std::vector<double> xCenters (Nx, 0.0);
  std::vector<double> yCenters (Ny, 0.0);

  for (int i = 0; i < Nx+1; i++){ xCoords[i] = i * dx; }
  for (int j = 0; j < Ny+1; j++){ yCoords[j] = j * dy; }
  for (int i = 0; i < Nx; i++){ xCenters[i] = 0.5 * (xCoords[i] + xCoords[i+1]); }
  for (int j = 0; j < Ny; j++){ yCenters[j] = 0.5 * (yCoords[j] + yCoords[j+1]); }

  struct BioParams {
    double gamma1 = 2.0;
    double gamma2 = 2.0;
    double tau_h = 0.7;   // maturation velocity scale
    double c1 = 0.68;
    double c2 = 0.08;
    double u_bar = 0.02;  // basal FSH scale
    double Lambda_bar = 0.1;
    double gamma_bar = 0.01; // source term spatial scale
    double Ys = 0.3;          // maturity threshold
  };

  struct HormoneParams {
    double U_min = 0.5;
    double c = 2.0;
    double M_ref = 4.5;
    double b1 = 0.08;
    double b2 = 2.25;
    double b3 = 1450.0;
  };

  // Variables
  int  totalFollicles = 1;
  std::vector<GridField2D> phi(totalFollicles);
  std::vector<GridField2D> g(totalFollicles);
  std::vector<GridField2D> h(totalFollicles);
  std::vector<GridField2D> Lambda(totalFollicles);
  std::vector<double> u_f(totalFollicles, 0.0);
  std::vector<double> m_f(totalFollicles, 0.0);
  double U = 0.0;
  double M = 0.0;
  std::vector<BioParams> bioParams(totalFollicles);
  for (int f=0; f < totalFollicles; ++f) {
    bioParams[f].Ys = Ys;
  }
  HormoneParams hormoneParams;

  // Resize fields;
  for (int f = 0; f < totalFollicles; ++f) {
    phi[f].resize(Nx, std::vector<double>(Ny, 0.0));
    g[f].resize(Nx, std::vector<double>(Ny, 0.0));
    h[f].resize(Nx, std::vector<double>(Ny, 0.0));
    Lambda[f].resize(Nx, std::vector<double>(Ny, 0.0));
  }

  auto filtered_less = yCoords | std::views::filter([Ys](double value) {return value < Ys;});
  std::vector<double> yLessThanYs(filtered_less.begin(), filtered_less.end());
  auto filtered_greater = yCoords | std::views::filter([Ys](double value) {return value >= Ys;});
  std::vector<double> yGreaterThanEqualToYs(filtered_greater.begin(), filtered_greater.end());

  double cx = 0.3;
  double cy = 0.20;
  double sigma = std::sqrt(0.002);

  // 1. Initialize Phi
  for (int f = 0; f < totalFollicles; ++f){
    for (int i = 0; i < cellsPerHalfCycle; ++i) {
      double x = xCenters[i];
      for (int j = 0; j < int(yLessThanYs.size()); ++j) {
        double y = yCenters[j];
        double dx_ = x -cx;
        double dy_ = y -cy;
        phi[f][i][j] = (1/(2 * std::numbers::pi * sigma * sigma)) * std::exp(-(dx_* dx_ + dy_ * dy_)/(2 * sigma * sigma));
        //phi[f][i][j] = evaluateGaussianQuadrature2D(
        //    xCoords[i], xCoords[i+1],
        //    yCoords[j], yCoords[j+1],
        //    cx, cy, sigma
        //    );
      }
    }
  }



  // 2. Compute maturity per follicle
  for (int f = 0; f < totalFollicles; ++f) {
    m_f[f] = 0.0; // TODO: Initialize at each time iteration
    for (int i = 0; i < Nx; ++i){
      for (int j = 0; j < Ny; ++j) {
        m_f[f] += yCenters[j] * phi[f][i][j] * dx * dy;
      }
    }
  }

  // 3. Compute Ovary maturity 
  M = 0; // TODO: Initialize at each time iteration
  for (int f = 0; f < totalFollicles; ++f) {
    M  += m_f[f];
  }

  // 4. Compute Global FSH
  U = 0.0;
  double numerator = 1.0 - hormoneParams.U_min;
  double denominator = 1.0 + std::exp(hormoneParams.c * (M - hormoneParams.M_ref));
  U = hormoneParams.U_min + (numerator / denominator);


  // 5. Compute local FSH
  for (int f = 0; f < totalFollicles; ++f) {
    double val = hormoneParams.b1 +  (std::exp(hormoneParams.b2 * m_f[f]) / hormoneParams.b3);
    if (val > 1.0) val = 1.0;
    u_f[f]  =  val * U;
  }


  // 6. Compute g, h, and Lambda
  // Omega1
  for (int f = 0; f < totalFollicles; ++f) {
    for (int cycle = 0; cycle < Nc; ++cycle) {
      int start = (cellsPerHalfCycle * 2 ) * cycle;
      int stop = cycle * (2 * cellsPerHalfCycle) + cellsPerHalfCycle;
      for (int  i = start ;  i < stop ; i++){
        int x_index = i; //int(xCoords[i]/dx);
        for (int j=0; j < int(yLessThanYs.size()); ++j) {
          double y = yCenters[j];
          int y_index = j; //int(yCoords[j]/dy);
          //g[f][x_index][y_index] = bioParams[f].gamma1 * u_f[f] + bioParams[f].gamma2;
          //h[f][x_index][y_index] = bioParams[f].tau_h * (-y * y + (bioParams[f].c1 * y + bioParams[f].c2) * (1.0 - std::exp(-u_f[f] / bioParams[f].u_bar)));
          //Lambda[f][x_index][y_index] = bioParams[f].Lambda_bar * std::exp(-(pow(y - bioParams[f].Ys, 2) / bioParams[f].gamma_bar)) * (1 - U);
          g[f][x_index][y_index] = 0.5;
          h[f][x_index][y_index] = 0.0;
          Lambda[f][x_index][y_index] = 0.0;
        }
      }
    }
  }

  // Omega2
  for (int f = 0; f < totalFollicles; ++f) {
    for (int cycle = 0; cycle < Nc; ++cycle) {
      int start = (cellsPerHalfCycle * 2 ) * cycle + cellsPerHalfCycle ;
      int stop = cycle * (2 * cellsPerHalfCycle) + cellsPerHalfCycle * 2;
      for (int  i = start; i < stop; ++i){
        int x_index = i; //int(xCoords[i]/dx);
        for (int j=0; j < int(yLessThanYs.size()); ++j) {
          int y_index = j; //int(yCoords[j]/dy);
          g[f][x_index][y_index] = 1.0;
          h[f][x_index][y_index] = 0.0;
          Lambda[f][x_index][y_index] = 0.0;
        }
      }
    }
  }

  // Omega 3
  for (int f = 0; f < totalFollicles; ++f) {
    for (int i = 0;  i < Nx; ++i) {
      int x_index = i; //int(xCoords[i]/dx);
                       // yGreaterThanEqualToYs contains the boundary value 1.0. So, use size()-1.
      for (int j = 0; j < int(yGreaterThanEqualToYs.size())-1; ++j){
        int y_index = static_cast<int>(std::round(yGreaterThanEqualToYs[j]/dy));
        double y = yCenters[y_index];
        //g[f][x_index][y_index] = 1.0;
        //h[f][x_index][y_index] = bioParams[f].tau_h * (- y * y + (bioParams[f].c1 * y + bioParams[f].c2) * (1.0 - std::exp(-u_f[f] / bioParams[f].u_bar)));
        //Lambda[f][x_index][y_index] = bioParams[f].Lambda_bar * std::exp(-(pow(y - bioParams[f].Ys, 2) / bioParams[f].gamma_bar)) * (1 - U);
          g[f][x_index][y_index] = 0.0;
          h[f][x_index][y_index] = 0.0;
          Lambda[f][x_index][y_index] = 0.0;
      }
    }
  }


  std::cout << std::string(60, '*') << "\n";
  std::cout
    << std::setw(15) << "m_f" 
    << std::setw(15) << "M"
    << std::setw(15) << "u_f"
    << std::setw(15) << "U" << "\n";

  std::cout
    << std::setw(15) << m_f[0]
    << std::setw(15) << M
    << std::setw(15) << u_f[0]
    << std::setw(15) << U  << "\n";
  std::cout << std::string(60, '*') << "\n";

  // Find interface indicies

  std::vector<int> interfaceXIndexContinuity(Nc, 0);
  std::vector<int> interfaceXIndexDoubling(Nc, 0);
  const int interfaceYIndexDirichlet = static_cast<int>(std::round(Ys /dy));
  std::vector<int> cellXIndexDirichlet(Nc * cellsPerHalfCycle, 0);
  for (int cycle = 0; cycle < Nc; ++cycle) {
    // Original
    interfaceXIndexContinuity[cycle] = cycle * 2 * cellsPerHalfCycle + cellsPerHalfCycle;
    interfaceXIndexDoubling[cycle] = (cycle + 1) * 2 * cellsPerHalfCycle; 
    // Switched
    //interfaceXIndexDoubling[cycle] = cycle * 2 * cellsPerHalfCycle + cellsPerHalfCycle;
    //interfaceXIndexContinuity[cycle] = (cycle + 1) * 2 * cellsPerHalfCycle; 

    //interfaceXIndexContinuity[cycle] = cycle * 2 * cellsPerHalfCycle + cellsPerHalfCycle-1;
    //interfaceXIndexDoubling[cycle] = (cycle + 1) * 2 * cellsPerHalfCycle - 1; 
    for (int idx = 0; idx < cellsPerHalfCycle; idx++){
      cellXIndexDirichlet[cycle * cellsPerHalfCycle + idx] = cycle * 2 * cellsPerHalfCycle + cellsPerHalfCycle + idx;
    }
  }

  //printVector("xCoordinates", xCoords);
  //printVector("xCenters", xCenters);
  //printVector("yCoordinates", yCoords);
  //printVector("yCenters", yCenters);
  //printVector("yLessThanYs", yLessThanYs);
  //printVector("yGreaterThanEqualToYs", yGreaterThanEqualToYs);
  printVector("interfaceXIndexContinuity", interfaceXIndexContinuity);
  printVector("interfaceXIndexDoubling", interfaceXIndexDoubling);
  //printVector("cellXIndexDirichlet", cellXIndexDirichlet);
  //std::cout << "interfaceYIndexDirichlet " << interfaceYIndexDirichlet << "\n";

  //printMatrix("phi0", phi[0]);
  //printMatrix("g0", g[0]);
  //printMatrix("h0", h[0]);
  //printMatrix("Lambda0", Lambda[0]);

  // Time-Stepping
  double dt;
  double cfl = 0.4;
  double t = 0.0;
  double t_final = 0.41;

  int step_count = 0;
  std::filesystem::create_directory("2dplain_results");
  std::string dirname = "2dplain_results";
  bool TIMESTEP = true;
  double mass = 0.0;
  while (t < t_final) {
    //for (int i = 0; i < Nx; ++i) {
    //  for (int j = 0; j < Ny; ++j) {
    //    mass += phi[0][i][j] * dx * dy;
    //  }
    //}
    mass = findAbsoluteMaximumInMatrix(phi[0]);
    std::cout << "mass " << mass << "\n";
    for (int f=0; f < totalFollicles; ++f) {
      // 0. Write output to file  
      std::string filename = dirname + "/step_" + std::to_string(step_count)+ ".csv";
      write_to_file(filename, xCenters, yCenters, phi[f]);
      // 1. Calculate "dt"
      double max_vel_x = findAbsoluteMaximumInMatrix(g[f]); 
      double max_vel_y = findAbsoluteMaximumInMatrix(h[f]); 
      double max_source = findAbsoluteMaximumInMatrix(Lambda[f]);
      double dt_vel_x = (cfl * dx ) /max_vel_x;
      double dt_vel_y = (cfl * dy) / max_vel_y;
      double dt_vel = std::min(dt_vel_x, dt_vel_y);
      // If righthand size is not equal to zero
      double dt_source = 0.0;
      if (max_source != 0.0) {
        dt_source = 1.0/max_source;
        dt = std::min(std::abs(dt_vel), std::abs(dt_source));
      } else {
        dt = std::abs(dt_vel);
      }
      if(t + dt > t_final) dt = t_final - t;
      std::cout
        << "Step " << step_count
        << " dt " << dt 
        << " max_g " << max_vel_x
        << " max_h " << max_vel_y
        << " max_l " << max_source
        << " cfl " << cfl 
        << " dx " << dx
        << " dy " << dy
        << " dt_g " << dt_vel_x
        << " dt_h " << dt_vel_y
        << " dt_vel " << dt_vel
        << " dt_source " << dt_source
        << "\n";
      // TODO: Remove Priting
      std::cout << " Nx " << Nx 
        << " Ny " << Ny
        << " Lx " << Lx
        << " Ly " << Ly
        << " Dc " << Dc
        << " Nc " << Nc
        << " Ys " << Ys 
        << " dx " << dx
        << " dy " << dy
        << " j_s "<< j_s
        << "\n";
      //std::cout << std::string(60, '*') << "\n";
      std::cout
        << std::setw(15) << "m_f" 
        << std::setw(15) << "M"
        << std::setw(15) << "u_f"
        << std::setw(15) << "U" << "\n";

      std::cout
        << std::setw(15) << m_f[0]
        << std::setw(15) << M
        << std::setw(15) << u_f[0]
        << std::setw(15) << U  << "\n";
      //std::cout << std::string(60, '*') << "\n";

      // 2. Update the function phi using RK3.
      GridField2D rhs1(Nx, std::vector<double>(Ny, 0.0));
      GridField2D rhs2(Nx, std::vector<double>(Ny, 0.0));
      GridField2D rhs3(Nx, std::vector<double>(Ny, 0.0));
      GridField2D phi_star1(Nx, std::vector<double>(Ny, 0.0));
      GridField2D phi_star2(Nx, std::vector<double>(Ny, 0.0));

      rhs1 = compute_rhs(Nx, Ny,
          phi[f], g[f], h[f], Lambda[f],
          interfaceXIndexContinuity,
          interfaceXIndexDoubling,
          cellXIndexDirichlet,
          interfaceYIndexDirichlet,
          dx, dy, dt
          );
      //printMatrix("phi at t=0", phi[f]);
      for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j){
          phi_star1[i][j] = phi[f][i][j] - rhs1[i][j];
        }
      }
      //printMatrix("phi_star1 at dt", phi_star1);
      rhs2 = compute_rhs(Nx, Ny,
          phi_star1, g[f], h[f], Lambda[f],
          interfaceXIndexContinuity,
          interfaceXIndexDoubling,
          cellXIndexDirichlet,
          interfaceYIndexDirichlet,
          dx, dy, dt
          );
      for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j){
          phi_star2[i][j] = phi[f][i][j] - (1.0/4.0) * rhs1[i][j] - (1.0/4.0) * rhs2[i][j]; 
        }
      }
      //printMatrix("phi_star2 at dt", phi_star2);
      rhs3 = compute_rhs(Nx, Ny,
          phi_star2, g[f], h[f], Lambda[f],
          interfaceXIndexContinuity,
          interfaceXIndexDoubling,
          cellXIndexDirichlet,
          interfaceYIndexDirichlet,
          dx, dy, dt
          );
      for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j){
          phi[f][i][j] = phi[f][i][j] 
            - (1.0/6.0) * rhs1[i][j]
            - (1.0/6.0) * rhs2[i][j]
            - (2.0/3.0) * rhs3[i][j]; 
        }
      }
      //printMatrix("phi at dt", phi[f]);

      // 3. Update the time step
      t += dt;
      ++step_count;
    }
    // 4. Update m_0, m_f, M, u_f, U, g, h, Lambda
    // 4.1 update maturity per follicle
    for (int f = 0; f < totalFollicles; ++f) {
      m_f[f] = 0.0; // TODO: Initialize at each time iteration
      for (int i = 0; i < Nx; ++i){
        for (int j = 0; j < Ny; ++j) {
          m_f[f] += yCenters[j] * phi[f][i][j] * dx * dy;
        }
      }
    }
    // 4.2. Compute Ovary maturity 
    M = 0; // TODO: Initialize at each time iteration
    for (int f = 0; f < totalFollicles; ++f) {
      M  += m_f[f];
    }

    // 4.3. Compute Global FSH
    U = 0.0;
    numerator = 1.0 - hormoneParams.U_min;
    denominator = 1.0 + std::exp(hormoneParams.c * (M - hormoneParams.M_ref));
    U = hormoneParams.U_min + (numerator / denominator);


    // 4.4 Compute local FSH
    for (int f = 0; f < totalFollicles; ++f) {
      double val = hormoneParams.b1 + (std::exp(hormoneParams.b2 * m_f[f]) / hormoneParams.b3);
      if (val > 1.0) val = 1.0;
      u_f[f]  =  val * U;
    }

    // 4.5. Compute g, h, and Lambda
    // Omega1
    for (int f = 0; f < totalFollicles; ++f) {
      for (int cycle = 0; cycle < Nc; ++cycle) {
        int start = (cellsPerHalfCycle * 2 ) * cycle;
        int stop = cycle * (2 * cellsPerHalfCycle) + cellsPerHalfCycle;
        for (int  i = start ;  i < stop ; i++){
          int x_index = i; //int(xCoords[i]/dx);
          for (int j=0; j < int(yLessThanYs.size()); ++j) {
            double y = yCenters[j];
            int y_index = j; //int(yCoords[j]/dy);
          g[f][x_index][y_index] = 0.5;
          h[f][x_index][y_index] = 0.0;
          Lambda[f][x_index][y_index] = 0.0;
            //g[f][x_index][y_index] = bioParams[f].gamma1 * u_f[f] + bioParams[f].gamma2;
            //h[f][x_index][y_index] = bioParams[f].tau_h * (-y * y + (bioParams[f].c1 * y + bioParams[f].c2) * (1.0 - std::exp(-u_f[f] / bioParams[f].u_bar)));
            //Lambda[f][x_index][y_index] = bioParams[f].Lambda_bar * std::exp(-(pow(y - bioParams[f].Ys, 2) / bioParams[f].gamma_bar)) * (1 - U);
          }
        }
      }
    }

    // Omega2
    for (int f = 0; f < totalFollicles; ++f) {
      for (int cycle = 0; cycle < Nc; ++cycle) {
        int start = (cellsPerHalfCycle * 2 ) * cycle + cellsPerHalfCycle ;
        int stop = cycle * (2 * cellsPerHalfCycle) + cellsPerHalfCycle * 2;
        for (int  i = start; i < stop; ++i){
          int x_index = i; //int(xCoords[i]/dx);
          for (int j=0; j < int(yLessThanYs.size()); ++j) {
            int y_index = j; //int(yCoords[j]/dy);
            g[f][x_index][y_index] = 1.0;
            h[f][x_index][y_index] = 0.0;
            Lambda[f][x_index][y_index] = 0.0;
          }
        }
      }
    }

    // Omega 3
    for (int f = 0; f < totalFollicles; ++f) {
      for (int i = 0;  i < Nx; ++i) {
        //int(xCoords[i]/dx);
        int x_index = i;
        // yGreaterThanEqualToYs contains the boundary value 1.0. So, use size()-1.
        for (int j = 0; j < int(yGreaterThanEqualToYs.size())-1; ++j){
          int y_index = static_cast<int>(std::round(yGreaterThanEqualToYs[j]/dy));
          double y = yCenters[y_index];
          //g[f][x_index][y_index] = 1.0;
          //h[f][x_index][y_index] = bioParams[f].tau_h * (- y * y + (bioParams[f].c1 * y + bioParams[f].c2) * (1.0 - std::exp(-u_f[f] / bioParams[f].u_bar)));
          //Lambda[f][x_index][y_index] = bioParams[f].Lambda_bar * std::exp(-(pow(y - bioParams[f].Ys, 2) / bioParams[f].gamma_bar)) * (1 - U);
          g[f][x_index][y_index] = 0.0;
          h[f][x_index][y_index] = 0.0;
          Lambda[f][x_index][y_index] = 0.0;
        }
      }
    }
    //if (TIMESTEP == true) break;
  }
  std::cout << "End of the program reached" << "\n";
  return 0;
}

