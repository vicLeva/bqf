
add_executable(random_uints
	random_uints.cpp
)

add_executable(ecoli_example
	ecoli.cpp
)

target_link_libraries(random_uints
	RSQF
)

target_link_libraries(ecoli_example
	RSQF
)

target_compile_definitions(ecoli_example PUBLIC DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
