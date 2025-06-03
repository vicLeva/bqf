cmake_minimum_required(VERSION 3.19)

project(RSQF)

SET (CMAKE_CXX_STANDARD 17)
SET(CMAKE_BUILD_TYPE Release)

enable_testing()

# C++ options
find_package(ZLIB REQUIRED)
add_compile_options(-Wall -Ofast -g -fpermissive -march=native -mtune=native) # for x86 based host
#add_compile_options(-Wall -Ofast -g -fpermissive -mcpu=native -mtune=native) # for ARM based host

# move all the binaries and libraries into dedicated directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(bundled/FQFeeder)
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(examples)
add_subdirectory(unitary_tests)

target_include_directories(RSQF PRIVATE ${CMAKE_SOURCE_DIR}/bundled/FQFeeder/include ${ZLIB_INCLUDE_DIRS})
target_link_libraries(RSQF PRIVATE fqfeeder ZLIB::ZLIB)