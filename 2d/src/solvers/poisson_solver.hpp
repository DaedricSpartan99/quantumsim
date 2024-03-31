#pragma once

#include "abstract_solver.hpp"
#include <memory>
#include "island_mesh.hpp"
#include "stiffness.hpp"
#include "load.hpp"

namespace qsim2d {

  /*
   * This class solves problems in the form:
   *
   *  div[rho(x) * grad[u(x)]] = f(x) :  x is internal
   *
   * with dirichlet zero boundary conditions:
   *
   *  u(x) = 0 : x on boundary
   */

  class PoissonSolver : public AbstractSolver {
    public:

      PoissonSolver(std::weak_ptr<IslandMesh> mesh, const ScalarField& rho, const ScalarField& f);

      virtual cpx_vector solve() override;

    private:

      // store scalarfield rho: stiffness
      StiffnessComponent stiffness;

      // store scalarfield f: load
      LoadComponent load;
      
      // store mesh 
      std::weak_ptr<IslandMesh> mesh;

      // store result
      cpx_vector solution;
  };
} 
