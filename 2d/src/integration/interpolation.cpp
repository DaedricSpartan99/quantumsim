#include "interpolation.hpp"

using namespace qsim2d;

Interpolator::Interpolator(const std::vector<InterpPair>& pairs) 
  : pairs(pairs)
{
}

double Interpolator::integrate(const ScalarField& f) const {
  
  double result = 0;

  for (auto pair : pairs) {
    result += pair.weight * f(pair.location);
  }

  return result;
}
