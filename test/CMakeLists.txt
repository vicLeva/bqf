
add_executable(quicktests
	quicktests.cpp
)
add_executable(resize_tests
	resize_tests.cpp
)

target_link_libraries(quicktests
	RSQF
)
target_link_libraries(resize_tests
	RSQF
)

target_compile_definitions(quicktests PUBLIC TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
target_compile_definitions(resize_tests PUBLIC TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")

