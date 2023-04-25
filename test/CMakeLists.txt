
add_executable(quicktests
	test.cpp
)

target_link_libraries(quicktests
	CQF
)

target_compile_definitions(quicktests PUBLIC TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")

