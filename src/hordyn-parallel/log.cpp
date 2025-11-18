#include "log.hpp"
#include <filesystem> // for creating directory

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
