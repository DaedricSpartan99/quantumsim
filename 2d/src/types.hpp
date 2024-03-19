#pragma once

#include <vector>
#include <array>
#include <complex>
#include <functional>

namespace qsim2d {
  
  // TODO use eigen  
  typedef std::complex<double> complex;
  typedef std::vector<complex> cpx_vector;
  typedef std::vector<cpx_vector> cpx_matrix;
  typedef std::array<double, 2> vertex_t;
  typedef std::array<vertex_t, 2> matrix_t;
  typedef std::function<vertex_t(vertex_t)> transform_t;
  typedef unsigned int index_t;
  typedef std::array<index_t, 3> triangle_t;
  
  /*
   * Operator overloading on vertex
   */

  // scalar product
  double operator*(const vertex_t&, const vertex_t&);
  
  // scalar vector multiplication
  vertex_t operator*(const vertex_t&, double);
  vertex_t operator*(double, const vertex_t&);
  
  // sum of vectors
  vertex_t operator+(const vertex_t&, const vertex_t&);

  /*
   * Operator overloading on vertex transformation matrix
   */

  // matrix vector multiplication
  vertex_t operator*(const matrix_t&, const vertex_t&);

  /*
   * Operator overloading on complex vector
   */

  // scalar product
  complex operator*(const cpx_vector&, const cpx_vector&);
  
  // scalar vector multiplication
  cpx_vector operator*(const cpx_vector&, complex);
  cpx_vector operator*(complex, const cpx_vector&);
  
  // sum of vectors
  cpx_vector operator+(const cpx_vector&, const cpx_vector&);

  /*
   * Operator overloading on complex matrix
   */

  cpx_vector operator*(const cpx_matrix&, const cpx_vector&);
}

