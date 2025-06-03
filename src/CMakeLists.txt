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

target_include_directories(RSQF
	PUBLIC
	"${PROJECT_SOURCE_DIR}/src/headers/"
)

target_link_libraries(bqf
	RSQF
)

set_property(TARGET RSQF PROPERTY CXX_STANDARD 17)
