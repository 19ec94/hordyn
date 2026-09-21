add_executable(main main.cpp domain.cpp computations.cpp compute_rhs.cpp bio_params.cpp
	compute_dt.cpp utils.cpp log.cpp)

target_compile_options(main PRIVATE
	-Wall
	-Wextra
	-Wpedantic
	-Warray-bounds
)
