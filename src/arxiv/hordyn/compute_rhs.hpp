#pragma once
#include "domain.hpp" // for domain class
#include "compute_dt.hpp" // for Time struct
#include <vector>

std::vector<double> compute_rhs(const Domain& domain,
				const Time& time,
				const std::vector<double>& phi,
				const std::vector<double>& g,
				const std::vector<double>& h,
				const std::vector<double>& source
				) ;
