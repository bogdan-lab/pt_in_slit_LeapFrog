cmake_minimum_required(VERSION 3.19)

project(LeapFrog_Pt_in_slit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#WARNINGS
add_library(full_set_warnings INTERFACE)
set(WARNINGS_AS_ERRORS FALSE)
include(./cmake/CompilerWarnings.cmake)
set_project_warnings(full_set_warnings)


add_subdirectory(src)


