#pragma once

#include "types.hpp"
#include <memory>

namespace qsim2d {

  /*
   * Class grouping all information needed to compute matrix components of finite element like stiffness or mass
   */
  
  class MatrixComponent {

    public:

      virtual matrix generate_matrix() const = 0;
  };

  /*
   * Class grouping all information needed to compute matrix components of finite element like stiffness or mass
   */
  
  class VectorComponent {

    public:

      virtual vector generate_vector() const = 0;
  };
}
