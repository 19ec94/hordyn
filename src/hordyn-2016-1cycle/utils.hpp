#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <iostream>

std::string create_folder(std::string path, int process_rank);

struct Log {
    std::vector<uint32_t> step_count;
    std::vector<double> dt;
    std::vector<double> time;
    std::vector<double> phi_max;
    std::vector<double> mass;
    std::vector<double> mf;
    std::vector<double> M;
    std::vector<double> U;
    std::vector<double> uf;

    // Method for writing logs to a csv file
    void write_summary(const std::string& filename) const {
        std::ofstream file(filename);
        file << "step,dt,time,phi_max,mass,local_maturity,global_maturity,global_fsh,local_fsh\n";
        for (size_t i = 0; i < step_count.size(); ++i) {
            file << step_count[i] << ","
                << dt[i] << ","
                << time[i] << ","
                << phi_max[i] << ","
                << mass[i] << ","
                << mf[i] << ","
                << M[i] << ","
                << U[i] << ","
                << uf[i] << "\n";
        }
    }

};


template<typename T>
void write_vector(const std::string& filename,
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
void write_matrix(const std::string& filename,
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
