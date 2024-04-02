#include "types.hpp"
#include "mesh.hpp"
#include "stiffness.hpp"
#include "load.hpp"
#include <iostream>
#include <vector>
#include "debug.hpp"

using namespace qsim2d;

matrix norm_fd_stiffness(size_t N);
vector fd_load(const std::vector<vertex_t>& vertices, std::function<double(vertex_t)> rhs);
double is_symmetric(const matrix& A);

template<class T>
T rhs_first_deg(T x, T y) {
  return x + y + 3; 
}

template<class T>
T rhs_second_deg(T x, T y) {
  return x * x + y * y;
}


template<class T>
T rhs_third_deg(T x, T y) {
  return x * x * x + y * y * y;
}

std::vector<InterpPair> GAUSS_FOUR_POINTS = {
        InterpPair{0.068464377, {0.112701665, 0.100000000}},
        InterpPair{0.109543004, {0.112701665, 0.443649167}},
        InterpPair{0.068464377, {0.112701665, 0.787298334}},
        InterpPair{0.061728395, {0.500000000, 0.056350832}},
        InterpPair{0.098765432, {0.500000000, 0.250000000}},
        InterpPair{0.061728395, {0.500000000, 0.443649167}},
        InterpPair{0.008696116, {0.887298334, 0.012701665}},
        InterpPair{0.013913785, {0.887298334, 0.056350832}},
        InterpPair{0.008696116, {0.887298334, 0.100000000}}
};

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


void output_vector(const vector& diff, const vector& v);

int main() {
  
  /*
   * Build a square [0, 1] x [0, 1] 
   */
  
  // square division per side
  const int N = 4;
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

  //npdebug("Number of vertices: ", vertices.size())
  //npdebug("Number of triangles: ", triangles.size())

  /*
   *  Perform finite difference equality test on stiffness
   */

  // construct simple mesh
  std::shared_ptr<IslandMesh> mesh = std::make_shared<IslandMesh>(vertices, triangles);

  // Build stiffness matrix and compute difference
  DirichletStiffness comp(mesh);

  matrix A_fe = comp.generate_matrix();

  matrix A_fd = norm_fd_stiffness(N-2);    // -2: exclude boundary

  matrix A_diff = A_fe - A_fd;

  /*
   * Output difference properties
   */
  std::cout << "Stiffness analysis" << std::endl << std::endl;

  // Size check
  std::cout << "FE Size: " << A_fe.col(0).size() << ", FD Size: " << A_fd.col(0).size() << std::endl;

  // Is-symmetric
  std::cout << "Is A symmetric: " << (is_symmetric(A_fe) < 1e-15) << std::endl;

  // norm
  std::cout << "Norm of difference: " << A_diff.norm() << std::endl;
  std::cout << std::endl;

  // diagonal
  std::cout << "Diagonal analysis" << std::endl;
  output_vector(A_fd.diagonal(0), A_fe.diagonal(0));
  std::cout << std::endl;

  // subdiagonals
  std::cout << "Inferior sub-diagonal analysis" << std::endl;
  output_vector(A_fd.diagonal(-1), A_fe.diagonal(-1));
  std::cout << std::endl;

  // sub-block diagonal
  std::cout << "First inferior sub-block-diagonal analysis" << std::endl;
  output_vector(A_fd.diagonal(-(N+1)), A_fe.diagonal(-(N+1)));
  std::cout << std::endl;
  std::cout << std::endl;

  /*
   *  Perform finite difference equality test on load
   */ 
  std::vector<vertex_t> internal_vert = mesh->get_internal_mesh().all_vertices();

  auto rhs_first = [](vertex_t v) -> double {
        return rhs_first_deg(v[0], v[1]);
      };

  DirichletLoad fe_load_comp(mesh, rhs_first);
  
  // get finite element load
  vector fe_l = fe_load_comp.generate_vector();
  
  // get finite difference load
  vector fd_l = fd_load(internal_vert, rhs_first) * h * h;

  // output properties
  std::cout << "Load analysis, rhs first order" << std::endl << std::endl;

  std::cout << "Vector equality test" << std::endl;
  output_vector(fd_l, fe_l);
  std::cout << std::endl;
  
  auto rhs_second = [](vertex_t v) -> double {
        return rhs_second_deg(v[0], v[1]);
      };

  auto rhs_third = [](vertex_t v) -> double {
        return rhs_third_deg(v[0], v[1]);
      };
  
  auto high_order_int = std::make_shared<const Interpolator>(GAUSS_MANY_POINTS);
  auto original_order_int = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS);

  DirichletLoad fe_load_second_comp(mesh, 
      rhs_second, high_order_int
  );
  
  // get finite element load
  vector fe_l_second = fe_load_second_comp.generate_vector();
  
  // get finite difference load
  vector fd_l_second = fd_load(internal_vert, rhs_second) * h * h;

  std::cout << "Load analysis, rhs second order" << std::endl << std::endl;

  std::cout << "Vector equality test" << std::endl;
  output_vector(fd_l_second, fe_l_second);
  std::cout << std::endl;

  const double exact_third_order = 1. / 10;
  const double high_order_test = high_order_int->integrate([](vertex_t) -> double { return 1; });
  const double numerical_third_order = high_order_int->integrate(rhs_third);
  const double numerical_third_order_original = original_order_int->integrate(rhs_third);

  std::cout << "Integration test second order" << std::endl;
  std::cout << "Integration of unit function: " << high_order_test << std::endl;
  std::cout << "Integral high order numerical value: " <<
    numerical_third_order << std::endl;
  std::cout << "Integral exact value: " <<
    exact_third_order << std::endl;
  std::cout << "Difference: " << abs(exact_third_order - numerical_third_order) << std::endl;
  std::cout << "Difference with original interpolation: " << abs(exact_third_order - numerical_third_order_original) << std::endl;
  std::cout << std::endl;
}

void output_vector(const vector& w, const vector& v) {

  vector diff(v - w);

  std::cout << "Difference norm: " << diff.norm() << std::endl;
  std::cout << "First three values ... end: " << 
    v[0] << " " <<
    v[1] << " " <<
    v[2] << " ... " <<
    v[v.size()-1] << std::endl;
  
  for (index_t i = 0; i < v.size(); ++i) {
    if (abs(diff[i]) > 1e-15)
      std::cout << "Difference at index " << 
        i << ", fd = " << 
        w[i] << ", fe = " <<
        v[i] <<
        std::endl;
  }
}

double is_symmetric(const matrix& A) {
  return (A - A.transpose()).norm();
}

matrix norm_fd_stiffness(size_t N) {
  
  const size_t size = (N+1) * (N+1);
  matrix A = matrix::Zero(size, size);

  for (int i = 0; i < N+1; ++i) {
    for (int j = 0; j < N+1; ++j) {

      // diagonal
      A(i*(N+1)+j,i*(N+1)+j) = 4.;

      // subdiagonal
      if (j < N)
        A(i*(N+1)+j,i*(N+1)+j+1) = -1.;

      if (j > 0)
        A(i*(N+1)+j,i*(N+1)+j-1) = -1.;

      // sub-block diagonal
      if (i < N)
        A(i*(N+1)+j,(i+1)*(N+1)+j) = -1.;

      if (i > 0)
        A(i*(N+1)+j,(i-1)*(N+1)+j) = -1.;
    }
  }

  return -A;
}


vector fd_load(const std::vector<vertex_t>& vertices, std::function<double(vertex_t)> rhs) {
  
  vector v = vector::Zero(vertices.size());

  for (size_t i = 0; i < vertices.size(); ++i) {
    v[i] = rhs(vertices[i]);
  }

  return v;
}
