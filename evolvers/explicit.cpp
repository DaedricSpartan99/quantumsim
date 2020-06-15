#include "explicit.hpp"

#include "grid/qsystem1D.hpp"
#include "grid/qsystem2D.hpp"

using namespace qsim::grid;
using namespace qsim::evo;

wave_vector explicit_scheme::evolve(const qsystem1D& system, double dt) const {
    using namespace std::complex_literals;
    return (wave_t(1.0) - (1i * dt / system.hbar()) * system.H()) * system.psi();
}

wave_grid explicit_scheme::evolve(const qsystem2D& system, double dt) const {
    // TODO
    return system.psi();
}
