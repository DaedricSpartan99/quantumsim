#pragma once

#include "abstract_solver.hpp"
#include <memory>
#include "island_mesh.hpp"
#include "stiffness.hpp"

namespace qsim2d {

  class ScalarField;

  /*
   * This class solves problems in the form:
   *
   *  div[rho(x) * grad[u(x)]] = 0 :  x is internal
   *
   * with dirichlet zero boundary conditions:
   *
   *  u(x) = 0 : x on boundary
   */

  class LaplaceSolver : public AbstractSolver {
    public:

      LaplaceSolver(std::weak_ptr<IslandMesh> mesh, const ScalarField& rho);

      virtual cpx_vector solve() override;

    private:

      // store scalarfield rho: stiffness
      StiffnessComponent stiffness;
      
      // store mesh 
      std::weak_ptr<IslandMesh> mesh;

      // lapack specific find matric kernel
      static std::vector<cpx_vector> find_kernel(const cpx_matrix&);
  };
} 
