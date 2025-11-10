cmake_minimum_required(VERSION 3.16)

project(hordyn VERSION 1.0 LANGUAGES CXX)

# C++ Standard festlegen (z.B C++20)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Speichere alle Executables unter ./bin Pfad.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Ordner include als Include-Pfad hinzufügen
# target_include_directories(include)

add_subdirectory(src/hordyn)
add_subdirectory(src/2d)
