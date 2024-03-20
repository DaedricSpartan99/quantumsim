#include "poisson_solver.hpp"
#include <eigen3/Eigen/QR>

using namespace qsim2d;

PoissonSolver::PoissonSolver(
    std::weak_ptr<IslandMesh> mesh, 
    const ScalarField& rho,
    const ScalarField& f
    ) 
  :
    stiffness(rho, mesh.lock()->get_internal_mesh()),
    load(f, mesh.lock()->get_internal_mesh()),
    mesh(mesh)
{
}

cpx_vector PoissonSolver::solve() {

  // determine system of equations
  cpx_matrix A = stiffness.generate_matrix();
  cpx_vector L = load.generate_matrix().col(0);

  // solve linear system
  return A.colPivHouseholderQr().solve(L);
}

