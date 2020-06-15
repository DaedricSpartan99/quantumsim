#pragma once

#include "fwd.hpp"
#include <complex>
#include <memory>
//#include "math/composition.hpp"
#include "evolver.hpp"

#include <cmath>

#include "debug.hpp"

namespace qsim {
    
    template <typename ... Coords>
    class potential;

    /*
     * Most general description of quantum system
     */
    template <class WaveFunction, typename ... Coords>
    class qsystem {
    protected:

        // the wave function
        WaveFunction wave;

        // an evolver
        std::shared_ptr<evolver> m_evolver;

    private:
        // mass
        double m;

        // a potential, external management
        std::shared_ptr<potential<Coords...>> pot;

        // plank constant
        double plank;

    public:
        qsystem(double _m, 
                const WaveFunction& _wave,
                std::shared_ptr<potential<Coords...>> _V,
                std::shared_ptr<evolver> _evolver,
                double _hbar = 1.0
                )
            : wave(_wave), m_evolver(_evolver), m(_m), pot(_V),  plank(_hbar) {

            if (pot == nullptr) {
                throw std::invalid_argument("The potential cannot be null");
            }
                
            if (m < 0)
                m = -m;
            
            if (plank < 0)
                plank = -plank;
        }
        
        // very important the virtual destructor
        virtual ~qsystem() {}

        // access to mass
        double mass() const {
            return m;
        }
        
        // eventually change the behaviour
        virtual void set_mass(double _m) {
            if (_m == 0)
                throw std::invalid_argument("Passed a null mass");
            this->m = std::abs(_m);
        }

        virtual void set_hbar(double hb) {
            npdebug("Hbar set: ", hb)
            if (hb == 0)
                throw std::invalid_argument("Passed a null plank constant");
            plank = std::abs(hb);
        }

        inline double hbar() const {
            return plank;
        }

        // evolution in time
        /*void evolve(double dt) {
            if (m_evolver != nullptr) {
                wave = std::move(m_evolver->evolve(*this, dt));
            }

            post(dt);
        }*/
        virtual void evolve(double) = 0;

        void set_evolver(std::shared_ptr<evolver> evo) {
            m_evolver = evo;
        }

        void set_potential(std::shared_ptr<potential<Coords...>> _V) {
            pot = _V;
        }
        
        // just in case after the evolution, some constraints must be set
        //virtual void post(double dt) {}


        inline const WaveFunction& psi() const {
            return wave;
        }

        inline void replace_wave(const WaveFunction& other) {
            wave = other;
        }

        inline void replace_wave(WaveFunction&& other) {
            wave = std::move(other);
        }

        inline const potential<Coords...>& V() const {
            return const_cast<const potential<Coords...>&>(*pot);
        }
        
        // averanges
        virtual double energy() const = 0;
        
        // normalization
        virtual double norm() const = 0;

        inline void normalize() {
            wave /= sqrt(norm());
        }

    protected:
        
        std::shared_ptr<potential<Coords...>> potential_ptr() const {
            return pot;
        }
    };
}
