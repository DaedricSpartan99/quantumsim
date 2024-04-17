#pragma once

#include "evolver.hpp"

namespace qsim::evo {

    class explicit_scheme : public qsim::evolver {

    public:

        virtual grid::wave_vector evolve(const grid::qsystem1D& system, double dt) const override;
        virtual grid::wave_grid evolve(const grid::qsystem2D& system, double dt) const override;
    };
}
