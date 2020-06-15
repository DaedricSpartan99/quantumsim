#pragma once

#include <vector>

namespace qsim::math {

    template <class T, class Obj>
    class ptr_composition {

        // identity (type safety)
        T identity; 
        
        // multiple of the identity
        T gain;
        
        // storage const pointers
        std::vector<const Obj*> components;

    public:

        ptr_composition(std::vector<const Obj*>&& comp, T id = 0.0, T _gain = 1.0)
            : identity(id), gain(_gain), components(std::move(comp)) {}

        ptr_composition(const ptr_composition&) = default;

        /*
         * Apply a linear operation over all elements
         * Requirements of V:
         *      - V& V::operator+=(V)
         *      - V operator*(T, In)
         */
        template<class V>
        V operator*(const V& input) const
        {
            V out((identity / gain) * input);

            for (const Obj* obj : components)
                out += (*obj) * input;

            return out *= gain;
        }

        qsim::math::ptr_composition<T, Obj>& operator*=(T mult) {
            identity *= mult;
            gain *= mult; 
            return *this;
        }

        qsim::math::ptr_composition<T, Obj>& operator+=(T add) {
            identity += add;
            return *this;
        }

        qsim::math::ptr_composition<T, Obj>& operator-=(T add) {
            identity -= add;
            return *this;
        }
        
        // self-inversion operator
        qsim::math::ptr_composition<T, Obj>& operator~() {
            identity *= -1.0;
            gain *= -1.0; 
            return *this;
        }
    };
}

template<typename T, class Obj, typename std::enable_if<std::is_floating_point<T>::value>::type>
qsim::math::ptr_composition<T, Obj> operator+(qsim::math::composition<T, Obj> input, T add) {
    // add to the first element
    return input += add;
}

template<typename T, class Obj>
qsim::math::ptr_composition<T, Obj> operator+(T add, qsim::math::ptr_composition<T, Obj> input) {
    // add to the first element
    return input += add;
}

template<typename T, class Obj>
qsim::math::ptr_composition<T, Obj> operator+(qsim::math::ptr_composition<T, Obj> input, T add) {
    // add to the first element
    return input += add;
}

template<typename T, class Obj>
qsim::math::ptr_composition<T, Obj> operator-(qsim::math::ptr_composition<T, Obj> input, T add) {
    // add to the first element
    return input -= add;
}

template<typename T, class Obj>
qsim::math::ptr_composition<T, Obj> operator-(T add, qsim::math::ptr_composition<T, Obj> input) {
    // add to the first element
    return ~(input -= add);
}

template<typename T, class Obj>
qsim::math::ptr_composition<T, Obj> operator*(qsim::math::ptr_composition<T, Obj> input, T mult) {
    // add to the first element
    return input *= mult;
}

template<typename T, class Obj>
qsim::math::ptr_composition<T, Obj> operator*(T mult, qsim::math::ptr_composition<T, Obj> input) {
    // add to the first element
    return input *= mult;
}


