#pragma once

#include <vector>
#include <array>
#include <complex>
#include <functional>
#include <eigen3/Eigen/Dense>

namespace qsim2d {
  
  // complex numbers
  typedef std::complex<double> complex;
  
  // real vectors
  typedef Eigen::VectorXd vector;
  typedef Eigen::MatrixXd matrix;
  
  // complex vectors
  typedef Eigen::VectorXcd cpx_vector;
  typedef Eigen::MatrixXcd cpx_matrix;

  // vertex and space transform
  typedef Eigen::Vector2d vertex_t;
  typedef Eigen::Matrix2d matrix_t;

  // index
  typedef Eigen::Index index_t;

  // functional
  typedef std::function<vertex_t(vertex_t)> transform_t;

  // triangles
  typedef std::array<index_t, 3> triangle_t;

  // fields
  typedef std::function<double(const vertex_t&)> ScalarField;
  typedef std::function<complex(const vertex_t&)> CpxScalarField;
  
  /*
   * Constexpr presets
   */  
  template<class T>
  inline const ScalarField unit_function = [](vertex_t) -> T {
    return 1.;
  };

  /*
   * Boundary conditions specification
   */
  enum BC {
    DIRICHLET,
    NEUMANN
  };
}

