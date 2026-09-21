#include "compute_dt.hpp"
#include "utils.hpp"
#include <algorithm>
#include <stdexcept>

double compute_dt(const Domain& domain, Time& time,
				const std::vector<double>& g,
				const std::vector<double>& h,
				const std::vector<double>& source) {
		double max_g = find_abs_max_value(g); 
		double max_h = find_abs_max_value(h);
		double max_source = find_abs_max_value(source);
		const double cfl = time.cfl;
		const double dx = domain.dx;
		const double dy = domain.dy;

		std::vector<double> dt_candidates;

		if(!is_zero(max_g)) {
				dt_candidates.push_back(cfl * (dx / max_g));
				time.max_g = max_g;
				time.dt_g = dt_candidates.back();
		}
		if (!is_zero(max_h)) {
				dt_candidates.push_back(cfl * (dy / max_h));
				time.max_h = max_h;
				time.dt_h = dt_candidates.back();
		}
		if (!is_zero(max_source)) {
				dt_candidates.push_back(1.0 / max_source);
				time.max_source = max_source;
				time.dt_source = dt_candidates.back();
		}

		if (dt_candidates.empty()) {
				throw std::runtime_error("compute_dt failed. dt = 0.0");
		}
		return *std::min_element(dt_candidates.begin(), dt_candidates.end());
}
