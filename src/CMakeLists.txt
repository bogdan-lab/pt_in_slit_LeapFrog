add_library(leapfrog_tracer SHARED
                spline.cpp
                particle.cpp
                settings.cpp
		)

target_link_libraries(leapfrog_tracer PUBLIC full_set_warnings m)
target_include_directories(leapfrog_tracer PUBLIC ../include)


find_package(OpenMP REQUIRED)

add_executable(program
                main.cpp)

target_link_libraries(program PRIVATE leapfrog_tracer OpenMP::OpenMP_CXX)
