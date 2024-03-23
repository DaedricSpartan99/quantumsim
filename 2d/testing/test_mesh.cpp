#include "types.hpp"
#include "mesh.hpp"
#include "island_mesh.hpp"
#include <iostream>
#include <vector>
#include <matplot/matplot.h>
#include <fstream>
#include <string>

using namespace qsim2d;

void output_mesh(const Mesh&, const std::string&);


int main() {
  
  /*
   * Build a square [0, 1] x [0, 1] 
   */
  
  // square division per side
  const int N = 5;

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

  // construct simple mesh 
  Mesh simple_mesh(vertices, triangles);

  // output simple mesh
  output_mesh(simple_mesh, "simple");
  
  // construct island mesh
  IslandMesh island_mesh(vertices, triangles);

  // output other mesh
  output_mesh(island_mesh.get_internal_mesh(), "internal");
}


void output_mesh(const Mesh& mesh, const std::string& file_prefix) {

  std::ofstream vert_stream(file_prefix + "_vertex.dat");

  for (const vertex_t& node : mesh.all_vertices()) {
    vert_stream << node[0] << " " << node[1] << std::endl;
  }

  vert_stream.close();

  std::ofstream triang_stream(file_prefix + "_triangles.dat");

  for (const triangle_t& triangle : mesh.all_triangles()) {
    triang_stream << triangle[0] << " " 
                << triangle[1] << " "
                << triangle[2] << std::endl;
  }

  triang_stream.close();
}


