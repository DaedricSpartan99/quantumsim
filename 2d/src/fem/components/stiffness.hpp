#pragma once

#include "matrix_component.hpp"

namespace qsim2d {

  /*
   * Class grouping all information needed to solve stiffness parameter
   * of a Poisson-like equation
   */
  
  class StiffnessComponent : public MatrixComponent {
    
    public:

      using MatrixComponent::MatrixComponent;

      virtual cpx_matrix generate_matrix() const override;
      
    private:

      // define gradients
      static constexpr std::array<point_t, 3> gradients = {
        {-1, -1}, 
        { 1, 0},
        { 0, 1}
      };
  };
}
