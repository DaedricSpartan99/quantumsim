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
    mesh(mesh),
    stiffness(mesh, rho),
    load(mesh, f)
{
}

vector PoissonSolver::solve() {

  // TODO: OPTIMIZE!!!!

  // determine system of equations
  npdebug("Generating stiffness")
  matrix A = stiffness.generate_matrix();
  npdebug("Stiffness matrix: ", A.diagonal())

  npdebug("Generating load")
  vector l = load.generate_vector();
  npdebug("Load vector: ", l)

  // solve linear system
  //solution = A.colPivHouseholderQr().solve(L);
  solution = A.fullPivLu().solve(l);

  // verify it's a solution
  npdebug("Relative error on solution: ", ((A * solution) - l).norm() / l.norm() * 100, "%")

  return solution;
}


