#include "types.hpp"
#include "mesh.hpp"
#include "island_mesh.hpp"
#include <iostream>
#include <vector>
#include <matplot/matplot.h>
#include <functional>
#include "solvers/poisson_solver.hpp"
#include <fstream>
#include "debug.hpp"

template<class T>
T exact(T x, T y) {
  return (1. - (x - 0.5)*(x - 0.5) / 4.) * (1. - (y - 0.5)*(y - 0.5) / 4.);
}

template<class T>
T rhs(T x, T y) {
  return -8. * (2. - 4. * ((x - 0.5) * (x - 0.5) + (y - 0.5)*(y - 0.5)) );
}

using namespace qsim2d;

void output_results(const cpx_vector& solution, const cpx_vector& exact);


int main() {
  
  /*
   * Build a square [0, 1] x [0, 1] 
   */
  
  // square division per side
  const int N = 10;
  const double h = 1. / N;

  // initialize components
  std::vector<vertex_t> vertices;
  std::vector<triangle_t> triangles;

  // draw first line of vertices 
  for (index_t j = 0; j < N+1; ++j)
    vertices.push_back(vertex_t{0., double(j) * h});
  
  // draw other vertices and bind triangles
  for (index_t i = 1; i < N+1; ++i) {

    // draw line of vertices 
    for (index_t j = 0; j < N+1; ++j)
      vertices.push_back(vertex_t{double(i) * h, double(j) * h});
    
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
  std::shared_ptr<IslandMesh> island_mesh = std::make_shared<IslandMesh>(vertices, triangles);

  // construct field
  ScalarField rhs_field([&](vertex_t v) -> complex {
        return (complex) rhs(v[0], v[1]);
      });

  // Build poisson solver
  npdebug("Is unit bound: ", static_cast<bool>(unit_function))
  npdebug("Is rhs bound: ", static_cast<bool>(rhs_field))
  PoissonSolver solver(island_mesh, unit_function, rhs_field);
  
  // solve system
  cpx_vector U = solver.solve();

  // compute exact solution for each vertex
  cpx_vector U_exact(U.size());
  const auto& internal_vert = island_mesh->get_internal_mesh().all_vertices();

  for (index_t i = 0; i < internal_vert.size(); ++i) {
    U_exact[i] = (complex) exact(vertices[i][0], vertices[i][1]);
  }

  // output to file
  output_results(U, U_exact);
}


void output_results(const cpx_vector& solution, const cpx_vector& exact) {
  
  cpx_vector sq_norm_diff = (solution - exact).cwiseAbs2();

  std::ofstream output("poisson.dat");

  output << std::setprecision(15);

  for (index_t i = 0; i < solution.size(); ++i) {
    output << solution[i].real() << " " << solution[i].imag() << " " << exact[i].real() << " " << exact[i].imag() << " " << sq_norm_diff[i].real() << std::endl;
  }

  output.close();
}
