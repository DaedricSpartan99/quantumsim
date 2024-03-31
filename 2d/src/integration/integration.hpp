#pragma once

#include "types.hpp"

namespace qsim2d {

  class Integrator {

    public:
      
      /*
       * Integrate over a fixed domain (specified in children classes)
       *
       * With the logic:
       *
       *  int{f(z)}
       */
      virtual double integrate(const ScalarField& f) const = 0;
  };
} 
