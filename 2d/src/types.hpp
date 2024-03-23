#pragma once

#include <vector>
#include <array>
#include <complex>
#include <functional>
#include <eigen3/Eigen/Dense>

namespace qsim2d {
  
  // TODO use eigen  
  typedef std::complex<double> complex;

  typedef Eigen::VectorXcd cpx_vector;
  typedef Eigen::MatrixXcd cpx_matrix;
  typedef Eigen::Vector2d vertex_t;
  typedef Eigen::Matrix2d matrix_t;
  typedef Eigen::Index index_t;

  typedef std::function<vertex_t(vertex_t)> transform_t;
  typedef std::array<index_t, 3> triangle_t;

  typedef std::function<complex(vertex_t)> ScalarField;
  
  /*
   * Constexpr presets
   */  
  const ScalarField unit_function;
}

