#pragma once

#include "evolver.hpp"
#include "grid/wave.hpp"

namespace qsim::evo {

    class crank_nicholson : public qsim::evolver {

    public:
        
        /*
         * Solve the half-implicit 1D scheme
         *
         * Complexity: O(N), N = grid size
         */
        virtual grid::wave_vector evolve(const grid::qsystem1D& system, double dt) const override;

        /*
         * Here is applied an approximation of the schroedinger equation to 
         * a (x,y) separable evolving equation, in order to
         * apply independent transformations
         *
         * Complexity: O(N*M), N = grid size along x, M = grid size along y
         */
        virtual grid::wave_grid evolve(const grid::qsystem2D& system, double dt) const override;
    };
}
