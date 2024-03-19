#pragma once

#include "types.hpp"

namespace qsim2d {
  
  class ScalarField {
    public:

      virtual complex evaluate(vertex_t) const = 0;

      inline complex operator()(vertex_t p) const {
        return evaluate(p);
      }
  };
}
