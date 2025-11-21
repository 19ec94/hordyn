#pragma once
#include <vector>
#include "domain.hpp"
#include "utils.hpp"
#include <algorithm>
#include <stdexcept>

struct Time {
		double dt = 0.0;
		double current_time = 0.0;
		double end_time = 0.0;
		double cfl = 0.0;
		uint32_t step_count = 0;

		// Info for post processing and cross-checking:
		double max_g = 0.0;
		double dt_g = 0.0;
		double max_h = 0.0;
		double dt_h = 0.0;
		double max_source = 0.0;
		double dt_source = 0.0;
};

double compute_dt(const Domain& domain, Time& time,
				const std::vector<double>& g,
				const std::vector<double>& h,
				const std::vector<double>& source);
