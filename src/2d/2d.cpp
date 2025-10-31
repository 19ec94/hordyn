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

using GridField2D = std::vector<std::vector<double>>;


// G1, SM, Differentiation phases
enum class CellPhase { Omega1, Omega2, Omega3};


std::ostream& operator<<(std::ostream& os, const CellPhase& phase) {
  switch (phase) {
    case CellPhase::Omega1: os << "G1"; break;
    case CellPhase::Omega2: os << "SM"; break;
    case CellPhase::Omega3: os << "DF"; break;
    default: os << "Unknown Phase"; break;
  }
  return os;
}


struct DomainParams {
  int Nc; // Number of cell cycles
  double Dc; // Duration of one cell cycle in x  (width in x)
  double Ys; // Maturity threshhold (0 < Ys < 1)

  int cellsPerHalfCycle; // Numer of cells per half-cycle in x

  int Nx; // Total number of cells in x (computed)
  int Ny; // Total number of cells in y (choosen s.t Ys aligns with Grid)

  double Lx; // Total domain length in x (computed)
  double Ly=1.0; // Total domain length in y (fixed)

  double dx; // Grid spacing in x (computed)
  double dy; // Grid spacing in y (computed)

  DomainParams(int Nc_, double Dc_, double Ys_, int cellsPerHalfCycle_, int Ny_)
    : Nc(Nc_), Dc(Dc_), Ys(Ys_), cellsPerHalfCycle(cellsPerHalfCycle_), Ny(Ny_)
  {
    Nx = 2 * Nc * cellsPerHalfCycle;
    Lx = Nc * Dc;
    dx = Lx / Nx;
    dy = Ly / Ny;

    double j_s = Ys / dy;
    if (std::abs(j_s - std::round(j_s)) > 1e-12) {
      throw std::runtime_error("Ys must align exactly with vertical grid lines");
    }
    if (Nx * Ny > 1e8) {
      throw std::runtime_error("Grid size too large for memory");
    }
  }
};


class Grid {
  public:
    const DomainParams& params;
    std::vector<double> xCoords; // edges in x, size Nx+1
    std::vector<double> yCoords; // edges in y, size Ny+1
    std::vector<double> xCenters; // centers in x, size Nx
    std::vector<double> yCenters; // centers in y, size Ny

    Grid(const DomainParams& p): params(p) {
      xCoords.resize(params.Nx + 1);
      yCoords.resize(params.Ny + 1);
      xCenters.resize(params.Nx);
      yCenters.resize(params.Ny);
      initialize();
    }

    void initialize() {
      for (int i=0; i <= params.Nx; ++i){
        xCoords[i] = i * params.dx;
      }
      for (int j=0; j <= params.Ny; ++j){
        yCoords[j] = j * params.dy;
      }
      for (int i=0; i < params.Nx; ++i){
        xCenters[i] = 0.5 * (xCoords[i] + xCoords[i+1]);
      }
      for (int j=0; j < params.Ny; ++j){
        yCenters[j] = 0.5 *(yCoords[j] + yCoords[j+1]);
      }
    }

    // Utility function: get cell index p for the horizontal cell cycle
    int cellCycleIndex(double x) const {
      int p = static_cast<int> (x / params.Dc);
      if (p < 0) p = 0;
      if (p >= params.Nc) p = params.Nc - 1;
      return p;
    }
};

class PhaseMap {
  public:
    const DomainParams& params;
    const Grid& grid;
    std::vector<std::vector<CellPhase>> phases; // size Nx x Ny

    PhaseMap(const DomainParams& p, const Grid& g) : params(p), grid(g) {
      phases.resize(params.Nx,
          std::vector<CellPhase> (params.Ny, CellPhase::Omega3));
      assignPhases();
    }

    void assignPhases() {
      for (int i = 0; i < params.Nx; ++i) {
        double x = grid.xCenters[i];
        int p = grid.cellCycleIndex(x);
        for (int j = 0; j < params.Ny; ++j){
          double y = grid.yCenters[j];
          if (y < params.Ys) {
            if (x >= p * params.Dc && x < p * params.Dc + 0.5 * params.Dc)
              phases[i][j] = CellPhase::Omega1;
            else if (x >= p * params.Dc + 0.5 * params.Dc && x < (p+1) * params.Dc)
              phases[i][j] = CellPhase::Omega2;
            else
              phases[i][j] = CellPhase::Omega3; // Outside subdomain ? Edge case
          }
          else {
            phases[i][j] = CellPhase::Omega3;
          }
        }
      }
    }

    CellPhase getPhase(int i, int j) const {
      return phases[i][j];
    }
};


enum class InterfaceType {
  None, 
  FluxContinuity, //Ω1 ↔ Ω2 interface with flux continuity
  FluxDoubling, //Ω2 ↔ Ω1 next cycle interface with flux doubling (mitosis) 
  Dirichlet //Ω2 ↔ Ω3 interface with Dirichlet condition (φ=0) 
};

std::ostream& operator<<(std::ostream& os, const InterfaceType& type) {
  switch (type) {
    case InterfaceType::None: os << "XX"; break;
    case InterfaceType::FluxContinuity: os << "FC"; break;
    case InterfaceType::FluxDoubling: os << "FD"; break;
    case InterfaceType::Dirichlet: os << "DR"; break;
    default: os << "Unknown Interface Type"; break;
  }
  return os;
}

struct Interface {
  int i; // x-index of interface (between cell i-1 and i)
  int j; // y-index of interface (between cell j-1 and j)

  InterfaceType type;
  //Aditional info (could be orientation, phase boundaries, etc)
  //For example, verticle or horizontal interface
  enum Orientation {Vertical, Horizontal } orientation;

  Interface(int i_, int j_, InterfaceType t, Orientation o) 
    : i(i_), j(j_), type(t), orientation(o) {}
};


struct InterfaceMaps {
  std::vector<std::vector<InterfaceType>> vertical;   // size: Nx+1 x Ny
  std::vector<std::vector<InterfaceType>> horizontal; // size: Nx x Ny+1

  InterfaceMaps(int Nx, int Ny) {
    vertical.resize(Nx + 1, std::vector<InterfaceType>(Ny, InterfaceType::None));
    horizontal.resize(Nx, std::vector<InterfaceType>(Ny + 1, InterfaceType::None));
  }
};


class InterfaceManager {
  const DomainParams& params;
  const Grid& grid;
  const PhaseMap& phaseMap;
  public:
  std::vector<Interface> verticalInterfaces;
  std::vector<Interface> horizontalInterfaces;

  InterfaceManager(const DomainParams& p, const Grid& g, const PhaseMap& ph)
    : params(p), grid(g), phaseMap(ph) {
      identifyInterfaces();
    }

  void identifyInterfaces() {
    // Identify vertical interfaces along x-edges between cells:
    // These exist at poisitons i where a change in phase or special
    // interface occurs.
    for (int i = 1; i < params.Nx; ++i){ //interface between cell i-1 and i
      for (int j = 0; j <params.Ny; ++j) {
        CellPhase leftPhase = phaseMap.getPhase(i-1, j);
        CellPhase rightPhase = phaseMap.getPhase(i,j);
        if (leftPhase != rightPhase) {
          // What to put here?
          InterfaceType type = InterfaceType::None;
          // Transmission condition: mitsosis end flux doubling (Omega2 to
          // Omega 1 of next cycle)
          if (leftPhase == CellPhase::Omega2 && rightPhase == CellPhase::Omega1) {
            // Transition end of cell cycle (mitosis)
            type = InterfaceType::FluxDoubling;
          } else if (leftPhase == CellPhase::Omega1 && rightPhase == CellPhase::Omega2) {
            // Transition start of cell cycle phase change
            type = InterfaceType::FluxContinuity;
          }
          // What other conditions can be added?
          if (type != InterfaceType::None) {
            verticalInterfaces.emplace_back(i, j, type, Interface::Vertical);
          }
        }
      }
    }
    // Special periodic boundary: last Omega2 (Nx-1) and first Omega1 (0)
    // Impose FluxDoubling condition
    for (int j = 0; j < params.Ny; ++j) {
      CellPhase leftPhase = phaseMap.getPhase(params.Nx - 1, j); // Last Omega2 column
      CellPhase rightPhase = phaseMap.getPhase(0, j); // first Omega1 column
      if (leftPhase != rightPhase) {
        InterfaceType type = InterfaceType::None;
        // Example: FluxDoubling between Omega2 -> Omega1 (end of cell cycle)
        if (leftPhase == CellPhase::Omega2 && rightPhase == CellPhase::Omega1) {
          type = InterfaceType::FluxDoubling;
        }
        else if (leftPhase == CellPhase::Omega1 && rightPhase == CellPhase::Omega2) {
          //type = InterfaceType::FluxContinuity; 
          type = InterfaceType::FluxDoubling; // because of periodicity
        }
        if (type != InterfaceType::None) {
          // note: use i = Nx and i=0 to represent periodic wrap-around 
          // (conceptually right of last Omega2, left of first Omega1)
          verticalInterfaces.emplace_back(params.Nx, j, type, Interface::Vertical);
          verticalInterfaces.emplace_back(0, j, type, Interface::Vertical);
        }
      }
    }


    // Identify horizontal interfaces where y = ys (maturity threshhold)
    // This will be at j = j_s (integer) where grid line matches Ys
    int j_s = static_cast<int>(std::round(params.Ys / params.dy));
    for (int i = 0; i < params.Nx; ++i) {
      // Interface between cells j_s-1 and j_s (horizontal interface)
      // This interface separates proliferate and differential areas
      CellPhase bottomPhase = phaseMap.getPhase(i, j_s-1);
      CellPhase topPhase = phaseMap.getPhase(i, j_s);
      if (bottomPhase == CellPhase::Omega2 && topPhase == CellPhase::Omega3) {
        horizontalInterfaces.emplace_back(i, j_s, InterfaceType::Dirichlet,
            Interface::Horizontal);
      }
    }
  }
};

InterfaceMaps makeInterfaceMaps(const DomainParams& params, const InterfaceManager& ifManager) {
  InterfaceMaps ifaceMaps(params.Nx, params.Ny);

  for (const auto& iface : ifManager.verticalInterfaces) {
    ifaceMaps.vertical[iface.i][iface.j] = iface.type;
  }
  for (const auto& iface : ifManager.horizontalInterfaces) {
    ifaceMaps.horizontal[iface.i][iface.j] = iface.type;
  }
  return ifaceMaps;
}

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


void initializePhi(const DomainParams& domainParams,
    const Grid& grid, const PhaseMap& phaseMap,
    std::vector<GridField2D>& phi, double cx, double cy, double sigma, int follicleIndex)
{
  // Set phi values to zero
  for (auto& row : phi[follicleIndex]) {
    std::fill(row.begin(), row.end(), 0.0);
  }

  for (int i = 0; i < domainParams.Nx; ++i) {
    double x = grid.xCenters[i];
    int p = static_cast<int>(x / domainParams.Dc);
    // Only initialize Omega1 cell in the first cell cycle (p == 0)
    if (p == 0) {
      for (int j = 0; j < domainParams.Ny; ++j) {
        double y = grid.yCenters[j];
        CellPhase phase = phaseMap.getPhase(i,j);
        if (phase == CellPhase::Omega1) {
          double dx = x - cx;
          double dy = y - cy;
          phi[follicleIndex][i][j] = std::exp(-(dx*dx + dy*dy) / (2 * sigma * sigma));
          //phi[follicleIndex][i][j] = evaluateGaussianQuadrature2D(
          //    grid.xCoords[i], grid.xCoords[i + 1],
          //    grid.yCoords[j], grid.yCoords[j + 1],
          //    cx, cy, sigma);
        }
      }
    }
  }
}


// Compute velocity fields g, h, and source term Lambda based on current phi, u_f, U
void computeGHLambda(const DomainParams& params, const Grid& grid, const PhaseMap& phaseMap,
    double u_f, double U,
    std::vector<GridField2D>& g,
    std::vector<GridField2D>& h,
    std::vector<GridField2D>& Lambda,
    const BioParams& bioParams, int follicleIndex) {

  for (int i = 0; i < params.Nx; ++i) {
    for (int j = 0; j < params.Ny; ++j) {
      double y = grid.yCenters[j];
      CellPhase phase = phaseMap.getPhase(i, j);
      switch (phase) {
        case CellPhase::Omega1:
          g[follicleIndex][i][j] = bioParams.gamma1 * u_f + bioParams.gamma2;
          h[follicleIndex][i][j] = bioParams.tau_h * ( - y * y + (bioParams.c1 * y + bioParams.c2) ) *
            (1.0 - std::exp(-u_f / bioParams.u_bar));
          Lambda[follicleIndex][i][j] = bioParams.Lambda_bar *
            std::exp(-(pow(y - bioParams.Ys, 2) / bioParams.gamma_bar)) * (1 - U);
          break;

        case CellPhase::Omega2:
          g[follicleIndex][i][j] = 1.0;
          h[follicleIndex][i][j] = 0.0;
          Lambda[follicleIndex][i][j] = 0.0;
          break;

        case CellPhase::Omega3:
          g[follicleIndex][i][j] = 1.0;
          h[follicleIndex][i][j] = bioParams.tau_h * (- y * y + (bioParams.c1 * y + bioParams.c2)) *
            (1.0 - std::exp(-u_f / bioParams.u_bar));
          Lambda[follicleIndex][i][j] = bioParams.Lambda_bar *
            std::exp(-(pow(y - bioParams.Ys, 2) / bioParams.gamma_bar)) * (1 - U);
          break;
      }
    }
  }
}

double computeFollicleMaturity(
    const DomainParams& domainParams,
    const Grid& grid,
    const GridField2D& phi)
{
  double result = 0.0; // Initialize sum to be zero before accumulation

  double cell_area = domainParams.dx * domainParams.dy;

  for (int i = 0; i < domainParams.Nx; ++i) {
    for (int j = 0; j < domainParams.Ny; ++j) {
      result += grid.yCenters[j] * phi[i][j] * cell_area;
    }
  }
  return result;
}

double computeOvaryMaturity(const std::vector<double>& m_f)
{
  double result = 0.0; // Initialize sum to be zero before accumulation
  for(auto mf_val : m_f) result += mf_val;
  return result;
}

double computeGlobalFSH(const double M, const HormoneParams& hormoneParams)
{
  double numerator = 1.0 - hormoneParams.U_min;
  double denominator = 1.0 + std::exp(hormoneParams.c * (M - hormoneParams.M_ref));
  return hormoneParams.U_min + ( numerator / denominator);
}

double computeLocalFSH(
    const double m_f,
    const double U,
    const HormoneParams& hormoneParams
    )
{
  double val = (hormoneParams.b1 + std::exp(hormoneParams.b2 * m_f)) / hormoneParams.b3;
  if(val > 1.0) val = 1.0;
  return val * U;
}


void writeGridDataToFile(const std::string& filename,
    const std::vector<double>& xCoords,
    const std::vector<double>& yCoords,
    const std::vector<double>& xCenters,
    const std::vector<double>& yCenters,
    const std::vector<std::vector<double>>& phi)
{
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file for writing\n";
    return;
  }

  // Write sizes on first line for convenience
  file << xCoords.size() << " " << yCoords.size() << " "
    << xCenters.size() << " " << yCenters.size() << " "
    << phi.size() << " " << (phi.empty() ? 0 : phi[0].size()) << "\n";

  // Write xCoords (Nx+1 values)
  for (auto val : xCoords)
    file << val << " ";
  file << "\n";

  // Write yCoords (Ny+1 values)
  for (auto val : yCoords)
    file << val << " ";
  file << "\n";

  // Write xCenters (Nx values)
  for (auto val : xCenters)
    file << val << " ";
  file << "\n";

  // Write yCenters (Ny values)
  for (auto val : yCenters)
    file << val << " ";
  file << "\n";

  // Write phi (Nx x Ny values)
  for (size_t i = 0; i < phi.size(); ++i) {
    for (size_t j = 0; j < phi[i].size(); ++j) {
      file << phi[i][j] << " ";
    }
    file << "\n";
  }

  file.close();
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

void updateStep( double dt,
    const GridField2D& g,
    const GridField2D& h,
    const GridField2D& Lambda,
    GridField2D& phi,
    const DomainParams& domainParams,
    const InterfaceMaps& ifaceMaps)
{
  int Nx = domainParams.Nx;
  int Ny = domainParams.Ny;
  GridField2D flux_left(Nx+1, std::vector<double>(Ny, 0.0));
  GridField2D flux_right(Nx+1, std::vector<double>(Ny, 0.0));
  GridField2D flux_y(Nx, std::vector<double>(Ny+1, 0.0));
}

int main() {
  DomainParams domainParams(8, 1.0, 0.75, 2, 8);

  Grid grid(domainParams);

  /*
     for (double value: grid.xCoords) {
     std::cout << value << " ";
     }
     std::cout << "\n";
     for (double value: grid.yCoords) {
     std::cout << value << " ";
     }
     std::cout << "\n";
     for (double value: grid.xCenters) {
     std::cout << value << " ";
     }
     std::cout << "\n";
     for (double value: grid.yCenters) {
     std::cout << value << " ";
     }
     std::cout << "\n";
     */


  PhaseMap phaseMap(domainParams, grid);

  /*
     std::cout << "Number of rows " << phaseMap.phases.size() << "\n";
     std::cout << "Number of columns " << phaseMap.phases[0].size() << "\n";

     for (int c = phaseMap.phases[0].size()-1; c>=0;  --c){
     for (int r = 0;  r < int(phaseMap.phases.size()); ++r){
     std::cout << phaseMap.phases[r][c] << " ";
     }
     std::cout << "\n";
     }
     */

  InterfaceManager interfaceManager(domainParams, grid, phaseMap);
  auto ifaceMaps = makeInterfaceMaps(domainParams, interfaceManager);

  // TODO: Remove Logging - print interface maps
  /*
     std::cout << "Vertical interfaces \n";
     for (int c = int(ifaceMaps.vertical[0].size())-1; c >=0;  --c) {
     for (int r = 0; r < int(ifaceMaps.vertical.size()); ++r) {
     std::cout << ifaceMaps.vertical[r][c] << " ";
     }
     std::cout << "\n";
     }
     std::cout << "Horizontal interfaces \n";
     for (int c = int(ifaceMaps.horizontal[0].size())-1; c >=0;  --c) {
     for (int r = 0; r < int(ifaceMaps.horizontal.size()); ++r) {
     std::cout << ifaceMaps.horizontal[r][c] << " ";
     }
     std::cout << "\n";
     }
     */

  /*
     for (const auto& iface : interfaceManager.verticalInterfaces) {
     double x_interface = grid.xCoords[iface.i];
     double y_center = grid.yCenters[iface.j];
     std::cout 
     << " at x = "
     << std::fixed << std::setprecision(3) << x_interface
     << ", y = "
     << std::fixed << std::setprecision(3) << y_center
     << " type = " << iface.type << "\n";
     }

     std::cout << "*******************************\n";

     for (const auto& iface : interfaceManager.horizontalInterfaces) {
     double y_interface = grid.yCoords[iface.j];
     double x_center = grid.xCenters[iface.i];
     std::cout
     << " at y = "
     << std::fixed << std::setprecision(3) << y_interface
     << ", x = "
     << std::fixed << std::setprecision(3) << x_center
     << " type = " << iface.type << "\n";
     }
     std::cout << "\n";
     */

  int numberOfFollicles = 1;
  std::vector<GridField2D> phi(numberOfFollicles);
  std::vector<GridField2D> g(numberOfFollicles);
  std::vector<GridField2D> h(numberOfFollicles);
  std::vector<GridField2D> Lambda(numberOfFollicles);
  std::vector<double> u_f(numberOfFollicles, 0.0);
  std::vector<double> m_f(numberOfFollicles, 0.0);
  // TODO: Make it time dependent vector or write to a file.
  double U = 0.0;
  double M = 0.0;
  std::vector<BioParams> bioParams(numberOfFollicles);
  for (int f = 0; f <numberOfFollicles; ++f) {
    bioParams[f].Ys = domainParams.Ys;
  }
  HormoneParams hormoneParams;

  // Resize fields for once
  for (int f = 0; f < numberOfFollicles; ++f) {
    phi[f].resize(domainParams.Nx, std::vector<double>(domainParams.Ny, 0.0));
    g[f].resize(domainParams.Nx, std::vector<double>(domainParams.Ny, 0.0));
    h[f].resize(domainParams.Nx, std::vector<double>(domainParams.Ny, 0.0));
    Lambda[f].resize(domainParams.Nx, std::vector<double>(domainParams.Ny, 0.0));
  }

  // 1. Initialize phi per follicle
  for (int f = 0; f < numberOfFollicles; ++f) {
    initializePhi(domainParams, grid, phaseMap, phi, 0.25, 0.15, std::sqrt(0.002), f); 
  }

  // 2. Calculate follicle maturity per follicle
  for (int f=0; f < numberOfFollicles; ++f) {
    m_f[f] = computeFollicleMaturity(domainParams, grid, phi[f]);
  }

  // 3. Calculate ovary maturity of all follicles
  M = computeOvaryMaturity(m_f);

  // 4. Compute global/plasma FSH 
  U = computeGlobalFSH(M, hormoneParams);

  // 5. Compute local FSH per follicle
  for (int f=0; f < numberOfFollicles; ++f) {
    u_f[f] = computeLocalFSH(m_f[f], U, hormoneParams);
  }

  // 6. Compute G, H, Lambda
  for (int f=0; f < numberOfFollicles; ++f) {
    computeGHLambda(domainParams, grid, phaseMap, u_f[f], U,
        g, h, Lambda, bioParams[f], f);
  }

  // TODO: Remove print Phi

  /*
     std::cout << "**************** Phi ******************\n";
      for (int c =  domainParams.Ny-1; c >=0; --c){
        for (int r = 0; r < domainParams.Nx; ++r) {
          std::cout << phi[0][r][c] << " ";
        }
        std::cout << "\n";
      }
     std::cout << "**************** g ******************\n";
     for (int c =  domainParams.Ny-1; c >=0; --c){
     for (int r = 0; r < domainParams.Nx; ++r) {
     std::cout << g[0][r][c] << " ";
     }
     std::cout << "\n";
     }
     std::cout << "****************** h *****************\n";

     for (int c =  domainParams.Ny-1; c >=0; --c){
     for (int r = 0; r < domainParams.Nx; ++r) {
     std::cout << h[0][r][c] << " ";
     }
     std::cout << "\n";
     }
     std::cout << "******************* Lambda ***************\n";

     for (int c =  domainParams.Ny-1; c >=0; --c){
     for (int r = 0; r < domainParams.Nx; ++r) {
     std::cout << Lambda[0][r][c] << " ";
     }
     std::cout << "\n";
     }
     */

  writeGridDataToFile("InitialData.txt", grid.xCoords, grid.yCoords, grid.xCenters,
      grid.yCenters, phi[0]);

  std::cout << std::string(60, '*') << "\n";

  std::cout
    << std::setw(15) << "m_f" 
    << std::setw(15) << "M"
    << std::setw(15) << "u_f"
    << std::setw(15) << "U" << "\n";

  std::cout << std::fixed << std::setprecision(10);

  std::cout
    << std::setw(15) << m_f[0]
    << std::setw(15) << M
    << std::setw(15) << u_f[0]
    << std::setw(15) << U  << "\n";
  std::cout << std::string(60, '*') << "\n";

  // Time-Stepping

  double dt;
  double cfl = 0.4;
  double t = 0.0;
  double t_final = 0.5;

  int step_count = 0;
  std::filesystem::create_directory("2d_full_model");

  bool TIMESTEP = true;
  while (TIMESTEP) {
    for (int f=0; f < numberOfFollicles; ++f){
      // 1. Calculate "dt"
      double max_vel_x = findAbsoluteMaximumInMatrix(g[f]); 
      double max_vel_y = findAbsoluteMaximumInMatrix(h[f]); 
      double max_source = findAbsoluteMaximumInMatrix(Lambda[f]);
      double dt_vel_x = (cfl * domainParams.dx ) /max_vel_x;
      double dt_vel_y = (cfl * domainParams.dy) / max_vel_y;
      double dt_vel = std::min(dt_vel_x, dt_vel_y);
      double dt_source = 1.0/max_source;
      dt = std::min(std::abs(dt_vel), std::abs(dt_source));
      if(t + dt > t_final) dt = t_final - t;

      // 2. Update the function phi 
      updateStep(dt, g[f], h[f], Lambda[f], phi[f], domainParams, ifaceMaps);
      // 3. Update the time step
      t += dt;
      ++step_count;
      // 4. Update m_0, m_f, M, u_f, U, g, h, Lambda
    }
    std::cout << "I stepped once" << "\n";
    TIMESTEP = false;
  }
  std::cout << std::string(60, '*') << "\n";
  std::cout << "End of program reached.\n";
  return 0;
}

