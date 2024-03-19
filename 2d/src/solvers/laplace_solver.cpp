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
  cpx_matrix solutions = LaplaceSolver::find_kernel(A);

  // take first solution
  return solutions.col(0);
}


cpx_matrix LaplaceSolver::find_kernel(const cpx_matrix& A) {

  using namespace Eigen;

  CompleteOrthogonalDecomposition<MatrixXcd> cod;
  cod.compute(A);

  // Find URV^T
  cpx_matrix V = cod.matrixZ().adjoint();
  cpx_matrix Null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
  cpx_matrix P = cod.colsPermutation();

  return P * Null_space; // Unpermute the columns
}
