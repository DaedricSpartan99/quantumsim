cmake_minimum_required(VERSION 3.0.0)
project (quantumsim2d)
#include_directories(grid evolvers)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-DNDEBUG)

include_directories(src fem/components fem/mesh)
add_library(quantumsim2d STATIC 
  src/types.cpp
  src/fem/components/mass.cpp
  src/fem/components/matrix_component.cpp
  src/fem/components/stiffness.cpp
  src/fem/mesh/mesh.cpp
  src/fem/mesh/island_mesh.cpp
  src/fem/mesh/poisson_solver.cpp
)

# vector/matrix test
add_executable(test_vector testing/test_vector.cpp)
target_link_libraries(test_vector quantumsim2d)
add_test(NAME test_vector
         COMMAND test_vector)

# mesh construction test
#add_executable(test_mesh testing/test_mesh.cpp)
#target_link_libraries(test_mesh quantumsim2d)
#add_test(NAME test_mesh
#         COMMAND test_mesh)

# stiffness matrix construction test
#add_executable(test_stiffness testing/test_stiffness.cpp)
#target_link_libraries(test_stiffness quantumsim2d)
#add_test(NAME test_stiffness
#         COMMAND test_stiffness)

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

# poisson dirichlet zero solver test
#add_executable(test_poisson testing/test_poisson.cpp)
#target_link_libraries(test_poisson quantumsim2d)
#add_test(NAME test_poisson
#         COMMAND test_poisson)

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


