find_package(MPI REQUIRED)

add_executable(main-2016 main.cpp domain.cpp utils.cpp time.cpp
	computations.cpp compute_rhs.cpp)

target_compile_options(main-2016 PRIVATE
	-Wall
	-Wextra
	-Wpedantic
	-Warray-bounds
)

target_link_libraries(main-2016 PRIVATE MPI::MPI_CXX)
