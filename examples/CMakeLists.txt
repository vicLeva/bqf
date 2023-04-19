
add_executable(random_uints
	random_uints.cpp
)

add_executable(user_example
	taro_oceans.cpp
)

target_link_libraries(random_uints
	CQF
)

target_link_libraries(user_example
	CQF
)
