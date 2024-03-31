#pragma once

#include <functional>
#include "matrix_component.hpp"

namespace qsim2d {

  /*
   * Class grouping all information needed to solve stiffness parameter
   * of a Poisson-like equation
   */
  
  class MassComponent : public MatrixComponent {
    
    public:

      using MatrixComponent::MatrixComponent;

      virtual cpx_matrix generate_matrix() const override;
      
    private:

      // Define basis functions
      static std::array<std::function<double(const vertex_t&)>, 3> basis;
  };
}
