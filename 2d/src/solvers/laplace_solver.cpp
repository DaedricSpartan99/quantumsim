#include "laplace_solver.hpp"
#include <eigen3/Eigen/QR>

using namespace qsim2d;

LaplaceSolver::LaplaceSolver(
    std::weak_ptr<IslandMesh> mesh, 
    const ScalarField& rho
    ) 
  :
    stiffness(rho, mesh.lock()->get_internal_mesh()),
    mesh(mesh)
{
}

cpx_vector LaplaceSolver::solve() {

  // determine system of equations
  cpx_matrix A = stiffness.generate_matrix();

  // find kernel of the matrix
  std::vector<cpx_vector> solutions = LaplaceSolver::find_kernel(A);

  // take first solution
  return solutions[0];
}


std::vector<cpx_vector> LaplaceSolver::find_kernel(const cpx_matrix& A) {

  using namespace Eigen;
  //typedef Matrix<qsim2d::complex, Dynamic, Dynamic> MatrixXcd;

  MatrixXcd eigen_A(A.size(), A[0].size());

  // copy from A 
  for (index_t j = 0; j < A[0].size(); ++j) {
    for (index_t i = 0; i < A.size(); ++i) {
      eigen_A(i,j) = A[i][j];
    }
  }

  CompleteOrthogonalDecomposition<MatrixXcd> cod;
  cod.compute(eigen_A);

  // Find URV^T
  MatrixXcd V = cod.matrixZ().adjoint();
  MatrixXcd Null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
  MatrixXcd P = cod.colsPermutation();
  Null_space = P * Null_space; // Unpermute the columns
  
  // convert to std::vector
  std::vector<cpx_vector> columns;

  for (int i = 0; i < Null_space.cols(); ++i) {
    
    // retrieve column
    VectorXcd col = Null_space.col(i);

    // copy data
    cpx_vector dest(col.size());
    VectorXcd::Map(&dest[0], col.size()) = col;

    // push into 
    columns.push_back(dest);
  }

  return columns;
}
