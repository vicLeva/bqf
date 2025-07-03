include(FetchContent)
FetchContent_Declare(gtest
  QUIET
  URL https://github.com/google/googletest/archive/refs/tags/v1.13.0.tar.gz
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gtest)

add_executable(bqf
	main.cpp
)

add_library(RSQF
	additional_methods.cpp
	rsqf.cpp
	abstract_bqf.cpp
	bqf_ec.cpp
	bqf_oom.cpp
	bqf_cf.cpp
)

target_link_libraries(RSQF PUBLIC gtest)

target_include_directories(RSQF
	PUBLIC
	"${PROJECT_SOURCE_DIR}/src/headers/"
)

target_link_libraries(bqf
	RSQF
)

set_property(TARGET RSQF PROPERTY CXX_STANDARD 17)
