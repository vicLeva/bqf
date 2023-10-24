add_library(RSQF
	additional_methods.cpp
	rsqf.cpp
	abstract_bqf.cpp
	bqf_ec.cpp
	bqf_oom.cpp
)

target_include_directories(RSQF
	PUBLIC
	"${PROJECT_SOURCE_DIR}/src/headers/"
)

set_property(TARGET RSQF PROPERTY CXX_STANDARD 17)
