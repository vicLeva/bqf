add_executable(resize_bench
	resize_bench.cpp
)

target_link_libraries(resize_bench
	RSQF
)

set_property(TARGET resize_bench PROPERTY CXX_STANDARD 17)
