#include "types.hpp"
#include "mesh.hpp"
#include "stiffness.hpp"
#include <iostream>
#include <vector>
#include "debug.hpp"

using namespace qsim2d;

Eigen::MatrixXcd fd_stiffness(size_t N);
double is_symmetric(const Eigen::MatrixXcd& A);

void output_vector(const Eigen::VectorXcd& diff, const Eigen::VectorXcd& v);

int main() {
  
  /*
   * Build a square [0, 1] x [0, 1] 
   */
  
  // square division per side
  const int N = 2;
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

  /*
   *  Perform finite difference equality test
   */

  // construct simple mesh
  Mesh mesh (vertices, triangles);

  // Build stiffness matrix and compute difference
  StiffnessComponent comp(unit_function, mesh);

  Eigen::MatrixXcd A_fe = comp.generate_matrix();
  Eigen::MatrixXcd A_fd = fd_stiffness(N);

  Eigen::MatrixXcd A_diff = A_fe - A_fd;

  /*
   * Output difference properties
   */

  // Is-symmetric
  std::cout << "Is A symmetric: " << is_symmetric(A_fe) << std::endl;

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

void output_vector(const Eigen::VectorXcd& w, const Eigen::VectorXcd& v) {

  Eigen::VectorXcd diff(v - w);

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

double is_symmetric(const Eigen::MatrixXcd& A) {
  return (A - A.transpose()).norm();
}

Eigen::MatrixXcd fd_stiffness(size_t N) {
  
  const size_t size = (N+1) * (N+1);
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(size, size);

  for (int i = 0; i < N+1; ++i) {
    for (int j = 0; j < N+1; ++j) {

      // diagonal
      A(i*N+j,i*N+j) = 4.;

      // subdiagonal
      if (j < N)
        A(i*N+j,i*N+j+1) = -1.;

      if (j > 0)
        A(i*N+j,i*N+j-1) = -1.;

      // sub-block diagonal
      if (i < N)
        A(i*N+j,(i+1)*N+j) = -1.;

      if (i > 0)
        A(i*N+j,(i-1)*N+j) = -1.;
    }
  }

  return A;
}