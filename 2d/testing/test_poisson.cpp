#include "types.hpp"
#include "mesh.hpp"
#include "island_mesh.hpp"
#include <iostream>
#include <vector>
#include <matplot/matplot.h>
#include <functional>

template<class T>
T exact(T x, T y) {
  return (1. - (x - 0.5)*(x - 0.5) / 4.) * (1. - (y - 0.5)*(y - 0.5) / 4.);
}

template<class T>
T rhs(T x, T y) {
  return -8. * (2. - 4. * ((x - 0.5) * (x - 0.5) + (y - 0.5)*(y - 0.5)) );
}

using namespace qsim2d;

int main() {
  
  /*
   * Build a square [0, 1] x [0, 1] 
   */
  
  // square division per side
  const int N = 20;

  // initialize components
  std::vector<vertex_t> vertices;
  std::vector<triangle_t> triangles;

  // draw first line of vertices 
  for (index_t j = 0; j < N+1; ++j)
    vertices.push_back(vertex_t{0., double(j) / N});
  
  // draw other vertices and bind triangles
  for (index_t i = 1; i < N+1; ++i) {

    // draw line of vertices 
    for (index_t j = 0; j < N+1; ++j)
      vertices.push_back(vertex_t{double(i) / N, double(j) / N});
    
    // construct two triangles per square (reference corner: bottom left)
    for (index_t j = 0; j < N; ++j) {

      index_t tl = (i-1) * (N+1) + j;
      index_t tr = (i-1) * (N+1) + j + 1;
      index_t bl = i * (N+1) + j;
      index_t br = i * (N+1) + j + 1;

      // bottom-left -> bottom-right -> top-right
      triangle_t T1 {bl, br, tr};

      // bottom-left -> top-left -> top-right
      triangle_t T2 {bl, tl, tr};

      // push them to collection
      triangles.push_back(T1);
      triangles.push_back(T2);
    }
  }

  using namespace matplot;

  // construct island mesh
  IslandMesh island_mesh(vertices, triangles);

  // Build poisson solver
}


