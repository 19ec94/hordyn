#pragma once

#include <fstream> // for writing output into file
#include <vector>
#include <string>

struct Log {
		std::vector<uint32_t> step_count_vec;
		std::vector<double> dt_vec;
		std::vector<double> time_vec;
		std::vector<double> phi_max_vec;
		std::vector<double> mass_vec;
		std::vector<double> local_maturity_vec;
		std::vector<double> global_maturity_vec;
		std::vector<double> global_fsh_vec;
		std::vector<double> local_fsh_vec;

		// Method for writing logs to a csv file
		void write_summary(const std::string& filename) const {
				std::ofstream file(filename);
				file << "step,dt,time,phi_max,mass,local_maturity,global_maturity,global_fsh,local_fsh\n";
				for (size_t i = 0; i < step_count_vec.size(); ++i) {
						file << step_count_vec[i] << ","
								<< dt_vec[i] << ","
								<< time_vec[i] << ","
								<< phi_max_vec[i] << ","
								<< mass_vec[i] << ","
								<< local_maturity_vec[i] << ","
								<< global_maturity_vec[i] << ","
								<< global_fsh_vec[i] << ","
								<< local_fsh_vec[i] << "\n";
				}
		}

		template<typename T>
				void write_vector_to_file(const std::string& filename,
								const std::vector<T>& vec) {
						std::ofstream file(filename);
						if (!file.is_open()) {
								throw std::runtime_error("Failed to open file: " + filename);
						}
						for (const auto& val : vec) {
								file << val << "\n";
						}
				}

		template<typename T>
				void write_matrix_to_file(const std::string& filename,
								const std::vector<T>& vec,
								size_t num_rows,
								size_t num_cols) {
						std::ofstream file(filename);
						if (!file.is_open()) {
								throw std::runtime_error("Cannot open file for writing: " + filename);
						}
						if (vec.size() != num_rows * num_cols) {
								throw std::runtime_error("Vector size does not match specified dimensions");
						}
						for (size_t row = 0; row < num_rows; ++row) {
								for (size_t col = 0; col < num_cols; ++col) {
										file << vec[row * num_cols + col];
										if (col + 1 < num_cols) file << " ";
								}
								file << "\n";
						}
				}
};



std::string create_folder(std::string path, int identifier);


