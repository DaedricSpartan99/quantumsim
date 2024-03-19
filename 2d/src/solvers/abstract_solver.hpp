#pragma once

#include <memory>
#include "types.hpp"

namespace qsim2d {
  
  class AbstractSolver {
    public:

      virtual cpx_vector solve() = 0;
  };

}
