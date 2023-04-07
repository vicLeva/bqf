add_library(CQF
	additional_methods.cpp
	filter.cpp
	backpack_cqf.cpp
)

target_include_directories(CQF
	PUBLIC
	"${PROJECT_SOURCE_DIR}/cqf/headers/"
)

set_property(TARGET CQF PROPERTY CXX_STANDARD 17)
