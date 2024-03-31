#pragma once

#include "integration.hpp"

namespace qsim2d {
  
  // Gauss pair (w, x)
  struct InterpPair {
    double weight;
    vertex_t location;    
  };

  const std::vector<InterpPair> GAUSS_TRI_POINTS = {
    InterpPair{1./6, {0, 0.5}},
    InterpPair{1./6, {0.5, 0}},
    InterpPair{1./6, {0.5, 0.5}}
  };

  class Interpolator : public Integrator {

    public:

      Interpolator(const std::vector<InterpPair>& pairs);

      virtual double integrate(const ScalarField& f) const override;

    private:

      std::vector<InterpPair> pairs;
  };
}
