cmake_minimum_required(VERSION 3.17)
project(Filters)

set(CMAKE_CXX_STANDARD 17)
find_package(iir)


add_executable(Filters main.cpp)
target_link_libraries(Filters iir::iir)
target_include_directories(Filters PUBLIC "/home/richard/CLionProjects/fft")