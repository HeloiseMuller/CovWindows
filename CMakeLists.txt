cmake_minimum_required(VERSION 3.9)
project(CovWindows LANGUAGES CXX)

find_package(OpenMP REQUIRED)

add_executable(CovWindows src/main.cpp)
target_compile_features(CovWindows PRIVATE cxx_std_17)
target_compile_options(CovWindows PRIVATE -W -Wall -Wextra -Wpedantic -Weffc++)
target_compile_options(CovWindows PRIVATE $<$<CONFIG:DEBUG>:-g>)
target_compile_options(CovWindows PRIVATE $<$<CONFIG:RELEASE>:-O2>)
target_link_libraries(CovWindows PRIVATE OpenMP::OpenMP_CXX)
