cmake_minimum_required(VERSION 3.19)

project(LeapFrog_Pt_in_slit
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#WARNINGS
add_library(full_set_warnings INTERFACE)
set(WARNINGS_AS_ERRORS FALSE)
include(./cmake/CompilerWarnings.cmake)
set_project_warnings(full_set_warnings)

#CONAN
find_program(CONAN conan)
if (NOT CONAN)
    message(FATAL_ERROR "Needs conan package manager: please visit https://conan.io and install")
endif ()
include(cmake/conan.cmake)
conan_check(VERSION 1.14.4 REQUIRED)
conan_cmake_run(REQUIRES
                catch2/2.13.3
        OPTIONS BASIC_SETUP CMAKE_TARGETS
        BUILD missing)

add_subdirectory(src)
add_subdirectory(tests)

