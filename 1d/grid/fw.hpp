#pragma once

#include <complex>

/*namespace qsim::math {
    template <typename T, class H>
    qsim::math::ptr_composition;
}*/

namespace qsim::grid {

    // grid coordinates
    typedef size_t grid_t;
    
    // hemiltonian class 
    //template<class H>
    //using grid_H = qsim::math::ptr_composition<double, H>;

    typedef std::complex<double> wave_t;

    class wave_vector;
}

