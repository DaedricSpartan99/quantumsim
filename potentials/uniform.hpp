#pragma once

#include "potential.hpp"

namespace qsim::pot {

    template <typename ...Coords> 
    class uniform : virtual public potential<Coords...> {

        double value; 

    public:

        uniform(double _value = 0) noexcept : value(_value) {}
        
        // set value
        inline void set_value(double val) {
            value = val;
        }

        // static cast to double
        inline double get_value() const {
            return value;
        }

        // static cast to double
        operator double() {
            return value;
        }
        
        virtual double operator()(const Coords&...) const override {
            return value;
        }
    };
}
