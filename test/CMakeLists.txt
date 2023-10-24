
add_executable(quicktests
	quicktests.cpp
)

target_link_libraries(quicktests
	RSQF
)

target_compile_definitions(quicktests PUBLIC TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")

