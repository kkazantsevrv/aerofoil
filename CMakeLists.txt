cmake_minimum_required(VERSION 3.0)

project(aerofoil)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(spline_lib STATIC spline.hpp spline.cpp utils.hpp)
add_executable(gidrofoil_test 
    #spline_test.cpp    
    gidrofoil_test.cpp
)

target_link_libraries(gidrofoil_test PRIVATE spline_lib)
