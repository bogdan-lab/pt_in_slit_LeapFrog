add_library(leapfrog_tracer SHARED
                spline.cpp
		)

target_link_libraries(leapfrog_tracer PUBLIC full_set_warnings m)
target_include_directories(leapfrog_tracer PUBLIC ../include)


add_executable(program
                main.cpp)

target_link_libraries(program PRIVATE leapfrog_tracer)
