#pragma once

#include <weak_ptr>
#include "types.hpp"

namespace qsim2d {
  
  class Mesh; 

  class FiniteElementSolver {
    public:

      FiniteElementSolver(std::weak_ptr<Mesh> mesh)
        : mesh(mesh) {}


      virtual cpx_vector solve() = 0;

    protected:

      std::weak_ptr<Mesh> mesh;
  };

}
