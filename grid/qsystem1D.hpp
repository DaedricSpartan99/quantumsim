#pragma once

#include "fwd.hpp"
#include "math/diagonals.hpp"
#include "math/diagonal.hpp"
#include "grid/wave.hpp"

#include "qsystem.hpp"
#include "math/composition.hpp"

namespace qsim::grid {

    typedef qsim::math::composition<
                              wave_t, 
                              math::diagonals<wave_t, 3>, 
                              math::diag_functor<wave_t>
                             > H_matrix_1D;

    typedef math::diagonals<wave_t, 2> P_matrix_1D;

    // concretization for a 1D grid
    class qsystem1D : public qgridsystem {
    private:

        // discretization
        double dx;

        // standard matrix A for the laplace operator in 1D
        static const math::diagonals<wave_t, 3> A;
        
        // centered differentiation application
        static const math::diagonals<wave_t, 2> P_templ;

        // determine first hamiltonian term as function of the mass and the discretization step
        qsim::math::diagonals<wave_t, 3> H_zero() const;

    public: 

        struct init_pack {
            std::function<qsim::wave_t (double)> f;
            size_t N;

            init_pack(const std::function<qsim::wave_t (double)>& _f = std::function<qsim::wave_t (double)>(),
                      size_t _N = 0
                     ) : f(_f), N(_N) {}

            inline qsim::wave_t operator*(double x) const {
                return f(x);
            }

            wave_vector generate(double dx) const;
        };

        qsystem1D(double _m, 
                  double dx,
                  std::shared_ptr<potential<size_t>> _V,
                  const init_pack& init = init_pack(),
                  std::shared_ptr<evolver> _evolver = nullptr,
                  double hbar = 1.0
                  );

        // implementations
        virtual double energy() const override;
        double position() const;
        double momentum() const;
        double probability(size_t, size_t) const;

        // normalize the wave function
        virtual double norm() const override;
        
        // operators 
        H_matrix_1D H() const;
        P_matrix_1D P() const;
        
        // evolution
        virtual void evolve(double dt) override;  
        
        void replace_wave(const init_pack&);

        inline size_t size() const {
            return wave.size();
        }

        /*
         * Discretization setter
         */

        void set_delta(double);
        double delta() const;

        /*
         * Access to boundaries
         * Idea: i = -1 -> x = 0
         */

        inline double x(size_t i) const {
            return dx * (i+1);
        }

        /*
         * iteration
         */

        wave_vector::iterator begin();
        wave_vector::iterator end();
        wave_vector::const_iterator begin() const;
        wave_vector::const_iterator end() const;
    };
}

