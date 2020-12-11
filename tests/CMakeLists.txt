enable_testing()

add_executable(SplineTests
                spline_tests.cpp)

target_link_libraries(SplineTests PRIVATE CONAN_PKG::catch2 leapfrog_tracer)


add_test(
    NAME spline_tests
    COMMAND SplineTests)
