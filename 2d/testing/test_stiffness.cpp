#include "types.hpp"
#include "mesh.hpp"
#include "stiffness.hpp"
#include <iostream>
#include <vector>
#include "debug.hpp"

using namespace qsim2d;

matrix norm_fd_stiffness(size_t N);
double is_symmetric(const matrix& A);

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
   *  Perform finite difference equality test
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

  
  // other parts (spoiler: should be zero)
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

  return A;
}
