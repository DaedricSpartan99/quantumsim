#pragma once

#include <cstddef>

#include <complex>

namespace qsim {

    /*
     * Define the WaveFunction internal type as std::complex<double>
     */
    
    typedef std::complex<double> wave_t;

    template <class WaveFunction, typename ...Coords>
    class qsystem;

    class evolver;

    namespace math {
        template <typename T>
        class matrix;
    }

    namespace grid {

        class wave_vector;
        typedef math::matrix<wave_t> wave_grid;

        class qsystem1D;
        class qsystem2D;

        typedef qsystem<wave_vector, size_t> qgridsystem;
        typedef qsystem<wave_grid, size_t, size_t> qbi_gridsystem;
    }
}
