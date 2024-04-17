cmake_minimum_required(VERSION 3.0.0)
project (quantumsim)
#include_directories(grid evolvers)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-DNDEBUG)

include_directories(. grid evolvers)
add_library(quantumsim STATIC grid/qsystem1D.cpp grid/qsystem2D.cpp grid/wave.cpp evolvers/explicit.cpp evolvers/crank_nicholson.cpp)