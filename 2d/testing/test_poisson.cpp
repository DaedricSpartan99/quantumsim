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
T g(T x) {
  return 1. - 4. * (x - 0.5)*(x - 0.5);
}

template<class T>
T exact(T x, T y) {
  return g(x) * g(y);
}

template<class T>
T rhs(T x, T y) {
  return -8. * ( g(x) + g(y) );
}

using namespace qsim2d;

std::vector<InterpPair> GAUSS_MANY_POINTS = {
  InterpPair{0.6846439e-01, {0.1063508e+00, 0.1063508e+00}},
  InterpPair{0.8563571e-01, {0.4718246e+00, 0.8452624e-01}},
  InterpPair{0.3858025e-01, {0.8372983e+00, 0.6270166e-01}},
  InterpPair{0.8563571e-01, {0.8452624e-01, 0.4718246e+00}},
  InterpPair{0.9876543e-01, {0.3750000e+00, 0.3750000e+00}},
  InterpPair{0.3782109e-01, {0.6654738e+00, 0.2781754e+00}},
  InterpPair{0.3858025e-01, {0.6270166e-01, 0.8372983e+00}},
  InterpPair{0.3782109e-01, {0.2781754e+00, 0.6654738e+00}},
  InterpPair{0.8696116e-02, {0.4936492e+00, 0.4936492e+00}}
};


void output_convergence(const std::vector<std::tuple<int, double, double>>& convergence);
void output_results(std::shared_ptr<IslandMesh>, const vector& solution, const vector& exact);
vector analytic_rhs(const std::vector<vertex_t>& vertices, double h);

std::shared_ptr<IslandMesh> generate_mesh(int N);
vector compute_exact(const std::vector<vertex_t>& vertices);

int main() {
  
  /*
   * Build a square [0, 1] x [0, 1] 
   */
  
  // square division per side
  const int N = 20;
  const double h = 1. / N;

    // construct island mesh
  std::shared_ptr<IslandMesh> island_mesh = generate_mesh(N);

  // construct field
  ScalarField rhs_field([&](vertex_t v) -> double {
        return rhs(v[0], v[1]);
      });

  // a finer integrator
  auto integrator = std::make_shared<const Interpolator>(GAUSS_MANY_POINTS);

  // Build poisson solver
  npdebug("Is unit bound: ", static_cast<bool>(unit_function<double>))
  npdebug("Is rhs bound: ", static_cast<bool>(rhs_field))
  PoissonSolver solver(island_mesh, unit_function<double>, rhs_field);
  
  // solve system
  vector U = solver.solve();

  // compute exact solution for each vertex
  const auto& internal_vert = island_mesh->get_internal_mesh().all_vertices();

  vector U_exact = compute_exact(internal_vert);

  // output to file
  output_results(island_mesh, U, U_exact);

  // print norm distance
  std::cout << "Norm of difference: " << (U - U_exact).norm() << std::endl;

  // Deepen on rhs
  std::cout << "Difference of rhs: " << (solver.get_stiffness().generate_matrix() * U_exact - analytic_rhs(internal_vert, h)).norm() << std::endl;

  // Rhs load
  std::cout << "Load difference on finite difference and finite element: " << (solver.get_load().generate_vector() - analytic_rhs(internal_vert, h)).norm() << std::endl;

  // Print rhs load
  std::cout << "RHS load, fd" << std::endl;
  Eigen::MatrixXd D(internal_vert.size(), 2);
  D << solver.get_load().generate_vector(), analytic_rhs(internal_vert, h);
  std::cout << D << std::endl;   

  /*
   * Perform convergence test
   */

  std::cout << std::endl << "Convergence test" << std::endl << std::endl;

  const int n_begin = 8;
  const int n_end = 40;
  
  std::vector<int> N_iters(n_end - n_begin);
  std::generate(N_iters.begin(), N_iters.end(), [n = n_begin]() mutable { return n++; });

  std::vector<std::tuple<int, double, double>> convergence;
  
  // simulate for each n
  for (auto n : N_iters) {

    std::cout << std::endl << "n = " << n << std::endl << std::endl;
    
    const double h = 1. / n;

    // construct island mesh
    std::shared_ptr<IslandMesh> island_mesh_n = generate_mesh(n);

    // construct field
    ScalarField rhs_field_n([&](vertex_t v) -> double {
      return rhs(v[0], v[1]);
    });

    // Build poisson solver
    PoissonSolver solver_n(island_mesh_n, unit_function<double>, rhs_field_n);
  
    // solve system
    vector U_n = solver_n.solve();

    // compute exact solution for each vertex
    const auto& internal_vert_n = island_mesh_n->get_internal_mesh().all_vertices();

    vector U_exact_n = compute_exact(internal_vert_n);

    convergence.push_back({n, h, (U_n - U_exact_n).norm()});
  }

  // output result
  output_convergence(convergence);
}


void output_convergence(const std::vector<std::tuple<int, double, double>>& convergence) {
  
  std::ofstream output("poisson_convergence.dat");

  output << std::setprecision(15);

  for (index_t i = 0; i < convergence.size(); ++i) {
    output << std::get<0>(convergence[i]) << " " << std::get<1>(convergence[i]) << " " << std::get<2>(convergence[i]) << std::endl;
  }

  output.close();
}

std::shared_ptr<IslandMesh> generate_mesh(int N) {

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

  // construct island mesh
  return std::make_shared<IslandMesh>(vertices, triangles);
}

vector compute_exact(const std::vector<vertex_t>& internal_vert) {

  vector U_exact(internal_vert.size());

  for (index_t i = 0; i < internal_vert.size(); ++i) {
    U_exact[i] = exact<double>(internal_vert[i][0], internal_vert[i][1]);
  }
  
  return U_exact;
}

vector analytic_rhs(const std::vector<vertex_t>& vertices, double h) {

  vector v = vector::Zero(vertices.size());

  for (int i = 0; i < vertices.size(); ++i) {
    v[i] = rhs(vertices[i][0], vertices[i][1]) * h * h;
  }

  return v;
}


void output_results(std::shared_ptr<IslandMesh> mesh, const vector& solution, const vector& exact) {

  const std::vector<vertex_t> vertices = mesh->get_internal_mesh().all_vertices();
  
  vector sq_norm_diff = (solution - exact).cwiseAbs2();

  std::ofstream output("poisson.dat");

  output << std::setprecision(15);

  for (index_t i = 0; i < solution.size(); ++i) {
    output << 
      vertices[i][0] << " " <<
      vertices[i][1] << " " <<
      solution[i] << " " << 
      exact[i] << " " << 
      sq_norm_diff[i] << std::endl;
  }

  output.close();
}
