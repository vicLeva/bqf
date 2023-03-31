
add_executable(benchmark
	benchmark.cpp
)

target_link_libraries(benchmark
	CQF
)

target_compile_definitions(benchmark PUBLIC BENCHMARK_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")