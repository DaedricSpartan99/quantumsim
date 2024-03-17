#include "fem_solver.hpp"
#include "unique_ptr"

namespace qsim2d {

  /*
   * This class solves problems in the form:
   *
   * rho(x) * Pois(u(x)) = lambda(x) u(x)
   *
   * with dirichlet zero boundary conditions.
   */

  class PoissonSolver : public FiniteElementSolver {
    public:

      PoissonSolver(const ScalarField&, std::weak_ptr<Mesh>);

      virtual cpx_vector solve() override;

    protected:
      
      // internal points
      std::unique_ptr<Mesh> internal_mesh;

      // boundary points
      std::vector<index_t> bnd_points;

    private:

      void separe_external_points();
  };
} 
