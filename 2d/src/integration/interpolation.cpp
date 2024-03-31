#include "interpolation.hpp"

using namespace qsim2d;

Interpolator::Interpolator(const std::vector<InterpPair>& pairs) 
  : pairs(pairs)
{
}

double Interpolator::integrate(const ScalarField& f, transform_t T) const {
  
  double result = 0;

  for (auto pair : pairs) {
    result += pair.weight * f(T(pair.location));
  }

  return result;
}
