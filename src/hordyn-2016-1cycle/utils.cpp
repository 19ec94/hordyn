#include "utils.hpp"


std::string create_folder(std::string path, int process_rank) {
    std::string folder_name = path + std::to_string(process_rank) + "/";
    std::filesystem::path directory_path(folder_name);

    try {
        std::filesystem::create_directories(directory_path);
        return directory_path.string();  // Return complete path as string
    } catch (const std::filesystem::filesystem_error& e) {
        throw std::runtime_error(std::string("Directory creation failed. ") + e.what());
    }
}

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

double find_abs_min_value(const std::vector<double>& values) {
		if (values.empty()) {
				throw std::runtime_error("Empty input vector, no minimum found");
		}
		auto min_iter = std::min_element(
						values.begin(), values.end(),
						[](double a, double b) { return std::abs(a) < std::abs(b); }
						);
		return std::abs(*min_iter);
}


bool is_zero(double value) {
		return std::fabs(value) < 1e-12; 
}
