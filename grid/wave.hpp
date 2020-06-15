#pragma once

#include <complex>
#include <vector>
#include <iosfwd>
#include "fw.hpp"

namespace qsim::grid {

    /*
    * Discretized wave function
    */
    class wave_vector : private std::vector<wave_t> {
    public:
        
        // inherit constructors
        using std::vector<wave_t>::vector;
        //wave_vector(size_t N = 0, wave_t value = wave_t(0));
        wave_vector(const wave_vector&);
        wave_vector(wave_vector&&);

        wave_vector& operator=(const wave_vector&);
        //wave_vector& operator=(wave_vector&&) = default;

        ~wave_vector();

        // inherit other methods
        using std::vector<wave_t>::size;
        using std::vector<wave_t>::clear;
        using std::vector<wave_t>::begin;
        using std::vector<wave_t>::end;
        using std::vector<wave_t>::front;
        using std::vector<wave_t>::back;
        using std::vector<wave_t>::operator[];

        using std::vector<wave_t>::push_back;
        using std::vector<wave_t>::insert;

        using iterator = std::vector<wave_t>::iterator;
        using const_iterator = std::vector<wave_t>::const_iterator;
        
        // define an out-of-bounds different behaviour for the const case
        /*inline wave_t operator[](int i) const {
            return at(i);
        }*/

        wave_t at(int i) const;

        wave_vector& operator+=(const wave_vector&);

        wave_vector& operator-=(const wave_vector&);
        wave_vector& operator*=(const wave_t&);
        wave_vector& operator/=(const wave_t&);

        /*
         * norm
         */

        double square_norm() const;

        /*
         * Add precompiled iterable object
         */

        template<class Iter>
        void push(const Iter& iter) {
            for (const auto& val : iter)
                this->push_back(val);
        }

        template<class Iterator>
        void push(Iterator begin, const Iterator end) {
            for (;begin != end; ++begin)
                this->push_back(*begin);
        }
    };
}

namespace qsim::grid {

wave_vector operator+(qsim::grid::wave_vector, const qsim::grid::wave_vector&);

qsim::grid::wave_vector operator-(qsim::grid::wave_vector, const qsim::grid::wave_vector&);

qsim::grid::wave_vector operator*(qsim::grid::wave_vector, const qsim::grid::wave_t&);
qsim::grid::wave_vector operator/(qsim::grid::wave_vector, const qsim::grid::wave_t&);
qsim::grid::wave_vector operator*(const qsim::grid::wave_t&, qsim::grid::wave_vector);
qsim::grid::wave_vector operator/(const qsim::grid::wave_t&, qsim::grid::wave_vector);

qsim::grid::wave_t operator*(qsim::grid::wave_vector, const qsim::grid::wave_vector&);

}

std::ostream& operator<<(std::ostream&, const qsim::grid::wave_vector&);

/*
 * Complex numbers operations
 */

namespace std {
    // complex conjugate
    qsim::grid::wave_vector conj(qsim::grid::wave_vector);
}

