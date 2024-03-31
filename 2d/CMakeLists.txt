cmake_minimum_required(VERSION 3.0.0)
project (quantumsim2d)
#include_directories(grid evolvers)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#add_definitions(-DNDEBUG)

include_directories(src src/fe src/fe/components src/fe/mesh src/integration src/solvers)
add_library(quantumsim2d STATIC 
  #src/fe/components/mass.cpp TODO
  src/fe/components/load.cpp
  src/fe/components/stiffness.cpp
  src/integration/interpolation.cpp
  src/fe/mesh/abstract_mesh.cpp
  src/fe/mesh/mesh.cpp
  src/fe/mesh/island_mesh.cpp
  #src/solvers/laplace_solver.cpp TODO
  src/solvers/poisson_solver.cpp 
)

# copy mesh tester into build directory
configure_file(${CMAKE_SOURCE_DIR}/testing/plot_mesh.py ${CMAKE_CURRENT_BINARY_DIR}/plot_mesh.py COPYONLY)

# mesh construction test
add_executable(test_mesh testing/test_mesh.cpp)
target_link_libraries(test_mesh quantumsim2d)
add_test(NAME test_mesh
         COMMAND test_mesh)

# stiffness matrix construction test
add_executable(test_stiffness testing/test_stiffness.cpp)
target_link_libraries(test_stiffness quantumsim2d)
add_test(NAME test_stiffness
         COMMAND test_stiffness)


# poisson dirichlet zero solver test
add_executable(test_poisson testing/test_poisson.cpp)
target_link_libraries(test_poisson quantumsim2d)
add_test(NAME test_poisson
         COMMAND test_poisson)


# mass matrix construction test
#add_executable(test_mass testing/test_mass.cpp)
#target_link_libraries(test_mass quantumsim2d)
#add_test(NAME test_mass
#         COMMAND test_mass)

# load vector construction test
#add_executable(test_load testing/test_load.cpp)
#target_link_libraries(test_load quantumsim2d)
#add_test(NAME test_load
#         COMMAND test_load)

# laplace dirichlet zero solver test
#add_executable(test_laplace testing/test_laplace.cpp)
#target_link_libraries(test_laplace quantumsim2d)
#add_test(NAME test_laplace
#         COMMAND test_laplace)

# heat dirichlet zero solver test
#add_executable(test_heat testing/test_heat.cpp)
#target_link_libraries(test_heat quantumsim2d)
#add_test(NAME test_heat
#         COMMAND test_heat)

# schrodinger dirichlet zero solver test
#add_executable(test_schrodinger testing/test_schrodinger.cpp)
#target_link_libraries(test_schrodinger quantumsim2d)
#add_test(NAME test_schrodinger
#         COMMAND test_schrodinger)


