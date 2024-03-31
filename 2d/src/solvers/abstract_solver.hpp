#pragma once

#include <memory>
#include "types.hpp"

namespace qsim2d {
  
  class AbstractSolver {
    public:
      
      // TODO: not enough general
      // TODO: uncomplete
      virtual cpx_vector solve() = 0;
  };

}
