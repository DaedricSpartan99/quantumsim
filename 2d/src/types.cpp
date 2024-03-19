#include "types.hpp"

using namespace qsim2d;

/*
 * Operator overloading on vertex
 */

// scalar product
double qsim2d::operator*(const vertex_t& a, const vertex_t& b) {

  return a[0] * b[0] + a[1] * a[1];
}

// scalar vector multiplication
vertex_t qsim2d::operator*(const vertex_t& v, double a) {
  
  return { v[0] * a, v[1] * a };
}

vertex_t qsim2d::operator*(double a, const vertex_t& v) {
  return v * a;
}

// sum of vectors
vertex_t qsim2d::operator+(const vertex_t& v, const vertex_t& w) {

  return { v[0] + w[0], v[1] + w[1] };
}

/*
 * Operator overloading on vertex transformation matrix
 */

// matrix vector multiplication
vertex_t qsim2d::operator*(const matrix_t& A, const vertex_t& x) {

  vertex_t b;

  for (int i = 0; i < 2; ++i) {

    b[i] = 0;

    for (int j = 0; j < 2; ++j) {
      b[i] += A[i][j] * x[j];
    }
  }

  return b;
}

/*
 * Operator overloading on complex vector
 */

// scalar product
complex qsim2d::operator*(const cpx_vector& v, const cpx_vector& w) {

  complex result = 0;
  const size_t N = std::min(v.size(), w.size());

  for (std::size_t i = 0; i < N; ++i) {
    result += v[i] * w[i];
  }

  return result;
}

// scalar vector multiplication
cpx_vector qsim2d::operator*(const cpx_vector& v, complex a) {

  cpx_vector result(v);

  for (std::size_t i = 0; i < v.size(); ++i) {
    result[i] *= a;
  }
  
  return result;
}

cpx_vector qsim2d::operator*(complex a, const cpx_vector& v) {
  return v * a;
}

// sum of vectors
cpx_vector qsim2d::operator+(const cpx_vector& v, const cpx_vector& w) {

  cpx_vector result(v);

  for (std::size_t i = 0; i < v.size(); ++i) {
    result[i] += w[i];
  }

  return result;
}

/*
 * Operator overloading on complex matrix
 */

cpx_vector qsim2d::operator*(const cpx_matrix& A, const cpx_vector& v) {

  cpx_vector b(A.size(), 0);

  for (int i = 0; i < A.size(); ++i) {

    for (int j = 0; j < v.size(); ++j)
      b[i] += A[i][j] * v[j];
  }

  return b;
}


