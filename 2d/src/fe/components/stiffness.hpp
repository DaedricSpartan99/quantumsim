#pragma once

#include "matrix_component.hpp"
#include "interpolation.hpp"
#include "fe_component.hpp"
#include "island_mesh.hpp"
#include "mesh.hpp"
#include <memory>

namespace qsim2d {

  /*
   * Class grouping all information needed to solve stiffness parameter
   * of a Poisson-like equation
   */
  
  class Stiffness : public MatrixComponent, public FEComponent {
    
    public:

      /*
       * Initialize the stiffness matrix with a custom field
       * and boundary conditions
       */
      Stiffness(
          std::weak_ptr<const AbstractMesh> mesh, 
          const ScalarField& field = unit_function<double>,
          std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS)
          );

      // implement matrix generation
      virtual matrix generate_matrix() const override;
 
    protected:

      // define gradients
      static const Eigen::Matrix<double, 2, 3> gradients;

      // integration
      std::shared_ptr<const Integrator> integrator;
  };

  /*
   * Specific to neumann boundary conditions
   */  
  class NeumannStiffness : public Stiffness {
    
    public:
      
      // Initialize the stiffness matrix with neumann boundary conditions
      NeumannStiffness(std::weak_ptr<const Mesh> mesh, 
          const ScalarField& field = unit_function<double>,
          std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS));

  };

  
  /*
   * Specific to dirichlet boundary conditions
   */  
  class DirichletStiffness : public Stiffness {
    
    public:
      
      // Initialize the stiffness matrix with dirichlet boundary conditions
      DirichletStiffness(std::weak_ptr<const IslandMesh> mesh, 
          const ScalarField& field = unit_function<double>,
          std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS));

  };

}
