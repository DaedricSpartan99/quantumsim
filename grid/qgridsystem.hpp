#pragma once

#include "qsystem.hpp"
#include "wave.hpp"
#include "math/ptr_composition.hpp"

namespace qsim::grid {

    // grid coordinates
    typedef size_t grid_t;
    
    // hemiltonian class 
    template<class H>
    using grid_H = qsim::math::ptr_composition<wave_t, H>;

    template<class H>
    class qgridsystem : public qsim::qsystem<size_t, wave_vector, grid_H<H>> {

    public:

        using qsystem<size_t, wave_vector, grid_H<H>>::qsystem;
        
        // do not copy real hemiltonian
        virtual grid_H<H> hemiltonian() const override {
            return grid_H<H>({hemiltonian_ptr()});
        }
        
        // this value depends on the exact grid
        virtual const H* hemiltonian_ptr() const = 0;

        inline size_t size() const {
            return this->psi().size();
        }
    };

}
