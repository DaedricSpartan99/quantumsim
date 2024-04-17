#pragma once

#include "matrix_component.hpp"
#include "interpolation.hpp"
#include "fe_component.hpp"
#include "island_mesh.hpp"
#include "mesh.hpp"
#include <memory>



namespace qsim2d {

  /*
   * Class grouping all information needed to solve mass parameter
   * of elliptic equations
   */
  
  class Mass : public MatrixComponent, public FEComponent {
    
    public:

      /*
       * Initialize the mass matrix with a custom field
       * and boundary conditions
       */
      Mass(
          std::weak_ptr<const AbstractMesh> mesh, 
          const ScalarField& field = unit_function<double>,
          std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS)
          );

      // implement matrix generation
      virtual matrix generate_matrix() const override;
 
    private:

      // Define basis functions
      static std::array<std::function<double(const vertex_t&)>, 3> basis;

      // integration
      std::shared_ptr<const Integrator> integrator;
  };

  /*
   * Specific to neumann boundary conditions
   */  
  class NeumannMass : public Mass {
    
    public:
      
      // Initialize the mass matrix with neumann boundary conditions
      NeumannMass(std::weak_ptr<const Mesh> mesh, 
          const ScalarField& field = unit_function<double>,
          std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS));

  };

  /*
   * Specific to dirichlet boundary conditions
   */  
  class DirichletMass : public Mass {
    
    public:
      
      // Initialize the mass matrix with dirichlet boundary conditions
      DirichletMass(std::weak_ptr<const IslandMesh> mesh, 
          const ScalarField& field = unit_function<double>,
          std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS));

  };

}
