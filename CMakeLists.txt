cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(msc)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(msc main.cpp)
add_executable(test_custom_struct test_custom_struct.cpp)
add_executable(test_1d_flat_vector test_1d_flat_vector.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

