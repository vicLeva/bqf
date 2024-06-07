
add_executable(quicktests
	quicktests.cpp
)
add_executable(resizeTests
	resizeTests.cpp
)
add_executable(perfResize
	perfResize.cpp
)

target_link_libraries(quicktests
	RSQF
)
target_link_libraries(resizeTests
	RSQF
)
target_link_libraries(perfResize
	RSQF
)

target_compile_definitions(quicktests PUBLIC TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
target_compile_definitions(resizeTests PUBLIC TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
target_compile_definitions(perfResize PUBLIC TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")

