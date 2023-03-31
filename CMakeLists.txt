cmake_minimum_required(VERSION 3.19)

project(CQF)

enable_testing()

# C++ options
add_compile_options(-Wall -Ofast -g -fpermissive) #-Werror

# move all the binaries and libraries into dedicated directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(cqf)
add_subdirectory(test)
add_subdirectory(examples)
add_subdirectory(unitary_tests)
add_subdirectory(benchmark)
