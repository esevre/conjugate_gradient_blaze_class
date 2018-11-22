cmake_minimum_required(VERSION 3.12)
project(cg_blaze_practice)

set(CMAKE_CXX_STANDARD 17)

find_package(blaze REQUIRED)


add_executable(cg_blaze_practice main.cpp matrix_setup.hpp conjugate_gradient.hpp)




