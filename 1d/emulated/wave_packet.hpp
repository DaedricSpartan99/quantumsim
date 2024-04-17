#pragma once

#include "region.hpp"

namespace qsim::emu {

    class observable;

    class wave_packet {
        
        // packet energy
        double _energy;
        
        // TODO abstract system
        const q_emusystem& system;

        std::shared_ptr<region> map; // from 0 to M
        std::shared_ptr<potential<size_t>> V; // from 0 to M

        std::vector<qsim::wave_t> A; // from 0 to M
        std::vector<qsim::wave_t> B; // from 0 to M
        
        // built-in observable integrals integrals
        static double norm(qsim::wave_t lambda, double x1, double x2, bool real);
        static double position(qsim::wave_t lambda, double x1, double x2, bool real);
        static double momentum(qsim::wave_t lambda, double x1, double x2, bool real);

    public:
        
        wave_packet(double energy, const q_emusystem& system);

        /*
         * iteration over regions
         */
        class iter_entry {
            const size_t k;
            const wave_packet& wave;
        public:

            iter_entry(size_t k, const wave_packet& wave);

            qsim::wave_t A() const;
            qsim::wave_t B() const;
            double V() const;

            // boundaries
            double x_min() const;
            double x_max() const;

            // check for real domain
            bool real() const;

            // obtain current domain lambda
            qsim::wave_t lambda() const;
        };

        class const_iterator {
            size_t k;
            const wave_packet& wave;
        public:
            const_iterator(const wave_packet&, size_t);

            const_iterator& operator++();
            const_iterator operator++(int);
            bool operator!=(const iterator&) const;

            iter_entry operator*() const;
        };

        const_iterator begin() const;
        const_iterator end() const;

        // regions
        inline size_t size() const {
            return map->size();
        }
    
        // lambda coefficient (complex)
        qsim::wave_t damp(size_t k) const;
        qsim::wave_t damp(double x) const;
        
        // packet energy setter
        void set_energy(double);
        double energy() const;

        /*
         * Access to A and B
         */
        qsim::wave_t coeff_A(size_t k) const;
        qsim::wave_t coeff_B(size_t k) const;

        qsim::wave_t coeff_A(double x) const;
        qsim::wave_t coeff_B(double x) const;

        inline std::pair<qsim::wave_t, qsim::wave_t> operator[](size_t k) const {
            return {coeff_A(k), coeff_B(k)};
        }
        
        /*
         * determine A and B in function of the potential
         * TODO, technical part
         */
        void eval_coefficients();

        /*
         * wave function access
         */
        qsim::wave_t psi(double x) const;

        // access to the boundary values
        // warning: psi(0) is undefined!!
        qsim::wave_t psi(size_t k) const;

        inline qsim::wave_t operator()(double x) const {
            return psi(x);
        }

        // wave function derivative, n = order
        qsim::wave_t derivative(double x, size_t n = 1) const;

        /*
         * integrals
         */
        // squared norm
        double norm() const;
        
        double position() const;
        double momentum() const;

        // custom observable action
        double observe(observable*) const;
    };
}
