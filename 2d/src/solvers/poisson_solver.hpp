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

  class PoissonSolver {
    public:

      PoissonSolver(std::weak_ptr<IslandMesh> mesh, const ScalarField& rho, const ScalarField& f, std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS));

      vector solve();

      const Stiffness& get_stiffness() const;

      const Load& get_load() const;

    private:

      // store mesh 
      std::weak_ptr<IslandMesh> mesh;

      // store scalarfield rho: stiffness
      Stiffness stiffness;

      // store scalarfield f: load
      Load load;

      // store result
      vector solution;
  };
} 
