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
      static constexpr std::array<std::function<double(const vertex_t&)>, 3> basis = {
        [](const vertex_t& z) -> double { return 1 - z[0] - z[1]; }, 
        [](const vertex_t& z) -> double { return z[0]; },
        [](const vertex_t& z) -> double { return z[1]; }
      };
  };
}
