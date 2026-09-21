#pragma once

#include <vector>
#include <iostream>

double find_abs_max_value(const std::vector<double>& values);
double find_abs_min_value(const std::vector<double>& values);
bool is_zero(double value);

template<typename T>
void print_matrix(const std::string& name, const std::vector<T>& matrix, int num_rows, int num_cols) {
		std::cout << "******** " << name << " *******"  << "\n";
		// Bottom row first
		for (int row = num_rows - 1; row >= 0; --row) {
				// Left to right
				for (int col = 0; col < num_cols; ++col) {
						std::cout << matrix[row * num_cols + col] << " ";
				}
				std::cout << "\n";
		}
}

template<typename T>
void print_vector(const std::string& name, const std::vector<T>& vector) {

		std::cout << "******** " << name << " *******"  << "\n";
		
		for (auto elem : vector) {
				std::cout << elem << " ";
		}
		std::cout << "\n";
}
