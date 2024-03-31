#pragma once

#include <functional>
#include "fe_component.hpp"
#include "interpolation.hpp"
#include "matrix_component.hpp"
#include "mesh.hpp"
#include "island_mesh.hpp"

namespace qsim2d {

  /*
   * Class grouping all information needed to solve stiffness parameter
   * of a Poisson-like equation
   */
  
  class Load : public VectorComponent, public FEComponent {
    
    public:

      /*
       * Initialize the stiffness matrix with a custom field
       * and boundary conditions
       */
      Load(
          std::weak_ptr<const AbstractMesh> mesh, 
          const ScalarField& field = unit_function<double>,
          std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS)
          );

      // implement matrix generation
      virtual vector generate_vector() const override;
 
    private:

      // Define basis functions
      static std::array<std::function<double(const vertex_t&)>, 3> basis;

      // integration
      std::shared_ptr<const Integrator> integrator;
  };

  /*
   * Specific to neumann boundary conditions
   */  
  class NeumannLoad : public Load {
    
    public:
      
      // Initialize the stiffness matrix with neumann boundary conditions
      NeumannLoad(std::weak_ptr<const Mesh> mesh, 
          const ScalarField& field = unit_function<double>,
          std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS));

  };

  
  /*
   * Specific to dirichlet boundary conditions
   */  
  class DirichletLoad : public Load {
    
    public:
      
      // Initialize the stiffness matrix with dirichlet boundary conditions
      DirichletLoad(std::weak_ptr<const IslandMesh> mesh, 
          const ScalarField& field = unit_function<double>,
          std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS));

  };

}
