cmake_minimum_required(VERSION 3.10)

option(BUILD_MACOS_INTEL "Build test programs" OFF)
option(BUILD_MACOS_ARM   "Build test programs" OFF)
option(BUILD_LINUX_INTEL "Build test programs" ON)
option(BUILD_LINUX_ARM   "Build test programs" OFF)

option(ENABLE_LTO    "Build test programs" OFF)
option(BUILD_PROFILE "Build test programs" OFF)
option(BUILD_DEBUG   "Build test programs" OFF)
option(BUILD_INFOS   "Build test programs" OFF)

#
#if(BUILD_MACOS_ARM)
#    message("BUILD_MACOS_ARM is not OFF/FALSE/0")
#endif()
#
#if(BUILD_MACOS_INTEL)
#    message("BUILD_MACOS_INTEL is not OFF/FALSE/0")
#endif()
#
#if(BUILD_MACOS_ARM OR BUILD_MACOS_INTEL)
#    if(BUILD_MACOS_ARM)
#        message("- BUILD_MACOS_ARM")
#    else()
#        message("- BUILD_MACOS_INTEL")
#    endif()
#endif()

SET(CMAKE_CXX_STANDARD     17)
SET(CMAKE_BUILD_TYPE  Release)

#
# MacOS with ARM M1/M2 processors
#

if(BUILD_MACOS_ARM)
    message("BUILD_MACOS AND BUILD_ARM => MACOS APPLE M1/M2")
    SET (CMAKE_C_COMPILER   /opt/homebrew/opt/llvm/bin/clang)     # LLVM 17
    SET (CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++)   # LLVM 17

    SET(CMAKE_CXX_FLAGS "-Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native -mtune=native")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math -fopenmp -std=c++17")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flax-vector-conversions") # gcc patch
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/opt/homebrew/Cellar/boost@1.76/1.76.0_5/include") # gcc patch

    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/lib")      # NemaTode
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
endif()

#
# MacOS with INTEL x86 processors
#

if(BUILD_MACOS_INTEL)
    message("BUILD_MACOS_INTEL configuration")
    SET(CMAKE_CXX_FLAGS "-Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_MKL_")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
endif ()

#
# Linux with INTEL x86 processors
#

if(BUILD_LINUX_INTEL)
    message("BUILD_LINUX_INTEL configuration")

#    SET(CMAKE_C_COMPILER   /usr/bin/clang-15)
#    SET(CMAKE_CXX_COMPILER /usr/bin/clang++-15)

    SET(CMAKE_CXX_FLAGS "-Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math -fopenmp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_MKL_")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g3")  #Uncomment to debug with valgrind
endif ()

#
# Linux with ARM processors
#

if(BUILD_LINUX_ARM)
    message("BUILD_LINUX_ARM configuration")
    #SET (CMAKE_C_COMPILER   /usr/bin/clang-15)
    #SET (CMAKE_CXX_COMPILER /usr/bin/clang++-15)
    SET (CMAKE_C_COMPILER   /opt/homebrew/bin/gcc-14)
    SET (CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-14)

    SET(CMAKE_CXX_FLAGS "-Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native -mtune=native")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g3 -ffast-math")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_MKL_")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")

#    SET(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -I/usr/local/include")
#    SET(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -I/opt/homebrew/include")
#    SET(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -I/opt/homebrew/Cellar/uhd/4.4.0.0_1/include/uhd")
#    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/lib")
#    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/opt/homebrew/lib")
#    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
endif ()

if(BUILD_PROFILE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass=loop-vectorize") #
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass-missed=loop-vectorize") #
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass-analysis=loop-vectorize") #
endif ()

if(ENABLE_LTO)
    message("LINK time optimization (LTO) is enable")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
endif ()


########################################################################################################################

project(bench_sorter)

#file(GLOB_RECURSE source_files ${CMAKE_CURRENT_SOURCE_DIR}/src/*)
file(GLOB_RECURSE source_files ./src/quicksort.cpp
        ./src/additional_methods.cpp
        src/sorting/std_4cores.cpp
        src/sorting/std_2cores.cpp
        ./src/tools/tools.cpp
        src/progress/progressbar.cpp
        ./src/bench_sorter.cpp)
add_executable(bench_sorter ${source_files})

target_link_libraries(bench_sorter pthread)
target_link_libraries(bench_sorter bz2)
target_link_libraries(bench_sorter curses)

########################################################################################################################

project(kmer_sorter)

#file(GLOB_RECURSE source_files ${CMAKE_CURRENT_SOURCE_DIR}/src/*)
file(GLOB_RECURSE source_files ./src/quicksort.cpp ./src/additional_methods.cpp
        ./src/progress/progressbar.cpp
        ./../src/additional_methods.cpp
        ./src/kmer_sorter.cpp)
add_executable(kmer_sorter ${source_files})

target_include_directories(kmer_sorter
        PUBLIC
        "${PROJECT_SOURCE_DIR}/../src/headers/"
)

target_link_libraries(kmer_sorter pthread)
target_link_libraries(kmer_sorter bz2)
target_link_libraries(kmer_sorter curses)

########################################################################################################################
