add_executable(2d 2d.cpp)
add_executable(2dplain 2dplain.cpp)

target_compile_options( 2d PRIVATE
	-Wall
	-Wextra
	-Wpedantic
	-Warray-bounds
)
target_compile_options( 2dplain PRIVATE
	-Wall
	-Wextra
	-Wpedantic
	-Warray-bounds
)
