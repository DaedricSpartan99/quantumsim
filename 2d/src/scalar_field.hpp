#pragma once

#include "types.hpp"

namespace qsim2d {
  
  class ScalarField {
    public:

      virtual std::complex evaluate(point_t) const = 0;

      inline std::complex operator()(point_t p) const {
        return evaluate(p);
      }
  };
}
