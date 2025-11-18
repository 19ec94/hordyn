find_package(MPI REQUIRED)

add_executable(main-parallel main.cpp domain.cpp computations.cpp compute_rhs.cpp 
	compute_dt.cpp utils.cpp log.cpp)

target_compile_options(main-parallel PRIVATE
	-Wall
	-Wextra
	-Wpedantic
	-Warray-bounds
)

target_link_libraries(main-parallel PRIVATE MPI::MPI_CXX)
