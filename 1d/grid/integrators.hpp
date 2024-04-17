#pragma once

#include "grid/wave.hpp"

namespace qsim::grid {
    
    template<class ... Coords>
    struct neighbourhood {

        static_assert((std::is_integral<Coords>::value && ...));

        const wave_vector& vector;
        std::tuple<Coords ...> jumps;
        
        neighbourhood(const wave_vector& _vector, Coords&& ... args)
            : vector(_vector), jumps(args ...) {}

        wave_t at(size_t offset, int k, const Coords& ... cs) const {
            if constexpr (sizeof...(Coords) == 0)
                return vector[static_cast<int>(offset) + k];
            else
                return std::apply([&](const Coords& ... Ns) { return vector[offset + k + ((Ns * cs) + ...) ]; }, jumps);
        }
    };
    
    // simplest integrator ever 
    template<class Op, class ... Jumps>    
    wave_t grid_integrate(const wave_vector& v, const Op& op, double dV, const Jumps& ... jumps) {
        wave_t out(0);
        neighbourhood<Jumps ...> map(v, std::forward(jumps) ...);

        for (size_t k = 0; k < v.size(); ++k)
            out += std::conj(v[k]) * op(map, k);

        return out *= dV;
    }
}
