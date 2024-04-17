#pragma once

#include "grid/wave.hpp"
#include "math/matrix.hpp"
#include "fwd.hpp"

namespace qsim {
    
    class evolver {
    public:

        virtual ~evolver() = default;
        
        virtual grid::wave_vector evolve(const grid::qsystem1D& system, double dt) const = 0;
        virtual grid::wave_grid evolve(const grid::qsystem2D& system, double dt) const = 0;
    };
}
