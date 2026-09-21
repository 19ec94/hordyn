#include "utils.hpp"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

// Important : it returns the absolute value of the maximum absolute value
double find_abs_max_value(const std::vector<double>& values) {
		if (values.empty()) {
				throw std::runtime_error("Empty input vector, no maximum found");
		}
		auto max_iter = std::max_element(
						values.begin(), values.end(),
						[](double a, double b) { return std::abs(a) < std::abs(b); }
						);
		return std::abs(*max_iter);
}


bool is_zero(double value) {
		return std::fabs(value) < 1e-12; 
}

