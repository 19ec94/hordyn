find_package(MPI REQUIRED)

add_executable(bridge-ode-pde main.cpp utils.cpp)

target_compile_options(bridge-ode-pde PRIVATE
	-Wall
	-Wextra
	-Wpedantic
	-Warray-bounds
)

target_link_libraries(bridge-ode-pde PRIVATE MPI::MPI_CXX)
# target_link_libraries(bridge-ode-pde PUBLIC matplot)
