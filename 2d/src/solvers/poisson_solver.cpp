#include "poisson_solver.hpp"
#include <eigen3/Eigen/QR>

#include "debug.hpp"

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
  npdebug("Generating stiffness")
  cpx_matrix A = stiffness.generate_matrix();
  npdebug("Stiffness matrix: ", A.diagonal().real())

  npdebug("Generating load")
  cpx_vector L = load.generate_matrix().col(0);
  npdebug("Load vector: ", L.real())

  // solve linear system
  //solution = A.colPivHouseholderQr().solve(L);
  solution = A.fullPivLu().solve(L);

  // verify it's a solution
  npdebug("Relative error on solution: ", ((A * solution) - L).norm() / L.norm() * 100, "%")

  return solution;
}


