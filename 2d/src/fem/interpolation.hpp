#pragma once

#include "types.hpp"

namespace qsim2d {
  
  // Gauss pair (w, x)
  struct GaussPair {
    double weight;
    vertex_t location;    
  };

  const std::vector<GaussPair> GAUSS_TRI_POINTS = {
    GaussPair{1./6, {0, 0.5}},
    GaussPair{1./6, {0.5, 0}},
    GaussPair{1./6, {0.5, 0.5}}
  };
}
