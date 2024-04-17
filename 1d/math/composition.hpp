#pragma once

#include <tuple>

#include "debug.hpp"

namespace qsim::math {
    
    /*
     * General sum composition of a mathematical object
     * This structure take advantages where properties
     * applies to all components singularly (ex: linear operations)
     *
     * T defines the field over the matrix is defined (integral, real, complex, ecc..)
     */
    template<typename T, class ...Obj>
    class composition {

        // identity (type safety)
        T identity; 
        
        // storage tuple
        std::tuple<Obj...> components;

    public:

        static constexpr size_t N = sizeof...(Obj);
        /*
         * id: first element, multiple of the identity
         * ...c: list of composing elements
         */
        composition(T id, const Obj& ...c) : identity(id), components(std::make_tuple(c...)) {
            //npdebug("Composition construction")
        }

        ~composition() {
            //npdebug("Destroy")
        }

        /*
         * Apply a linear operation over all elements
         * Requirements of V:
         *      - V operator+(V, V)
         *      - V operator*(T, V)
         */
        template<class V>
        V operator*(const V& input) const
        {
            V out(identity * input);

            if constexpr (sizeof...(Obj) > 0) {
                out += std::apply([&](const Obj& ...objs) -> V { return ((objs * input) + ...); }, components);
            }

            return out;
        }

        /*template<class V>
        inline V operator<<(const V& input) const {
            // TODO, do it better
            return composition::operator*(input);
        }*/

        /*
         * Extend the composition by pushing another element at the end
         */ 
        template<class Another>
        qsim::math::composition<T, Obj..., Another> extend(Another&& addon) const {
            return qsim::math::composition<T, Obj..., Another>(identity, (std::get<Obj>(components), ...), addon);
        }
       
        /*
         * Get element by positional index
         * index 0 = identity part
         */ 
        template <size_t K>
        auto& get() {
            if constexpr (K == 0)
                return identity;
            else {
                return std::get<K-1>(components);
            }
        }

        template <size_t K>
        const auto& get() const {
            if constexpr (K == 0)
                return identity;
            else {
                return std::get<K>(components);
            }
        }
        
        /*
         * Operator modifiers
         */ 
        qsim::math::composition<T, Obj...>& operator*=(T mult) {
            identity *= mult;
            if constexpr (sizeof...(Obj) > 0)
                std::apply([&mult](Obj& ...obj) noexcept { (..., (obj *= mult)); }, components);
            return *this;
        }

        qsim::math::composition<T, Obj...>& operator+=(T add) {
            identity += add;
            return *this;
        }

        qsim::math::composition<T, Obj...>& operator-=(T add) {
            identity -= add;
            return *this;
        }

        qsim::math::composition<T, Obj...>& operator~() {
            return (*this) *= T(-1);
        }

        /*
         * Random access
         */

        T at(size_t i, size_t j) const {
            T out((i == j) ? identity : T(0));
            return out += std::apply([&](const Obj& ...obj) noexcept { return (obj.at(i,j) + ...); }, components);
        }
    };
}

namespace qsim::math {

/*
 * double and complex addiction (specialization)
 */
template<typename T, class ...Obj>
qsim::math::composition<T, Obj...> operator+(qsim::math::composition<T, Obj...> input, T add) {
    // add to the first element
    return input += add;
}

template<typename T, class ...Obj>
qsim::math::composition<T, Obj...> operator+(T add, qsim::math::composition<T, Obj...> input) {
    // add to the first element
    return input += add;
}

template<typename T, class ...Obj>
qsim::math::composition<T, Obj...> operator-(qsim::math::composition<T, Obj...> input, T add) {
    // add to the first element
    return input -= add;
}

template<typename T, class ...Obj>
qsim::math::composition<T, Obj...> operator-(T add, qsim::math::composition<T, Obj...> input) {
    // add to the first element
    return ~(input -= add);
}

/*
 * By scalar multiplication
 */


template<typename T, class ...Obj>
qsim::math::composition<T, Obj...> operator*(qsim::math::composition<T, Obj...> input, T mult) {
    // add to the first element
    return input *= mult;
}

template<typename T, class ...Obj>
qsim::math::composition<T, Obj...> operator*(T mult, qsim::math::composition<T, Obj...> input) {
    // add to the first element
    return input *= mult;
}

/*
 * Push element into a new composition if Another is not a floating point
 */
template<typename T, class ...Obj, class Another, typename = typename std::enable_if<!std::is_same<T, Another>::value>::type>
qsim::math::composition<T, Obj..., Another> operator+(qsim::math::composition<T, Obj...>&& input, Another&& addon) {
    // expand to tuple and push the new element
    return qsim::math::composition<T, Obj..., Another>(input.identity, (std::get<Obj>(input.components), ...), addon);
}

}

/*
 * Element-wise multiplication
 */
/*namespace qsim::helper {
    template<class Op, class V>
    struct application {
        V operator()(const Op& op, const V& v) {
            return op * v;
        }
    };
}

template<class ...Obj, class V>
V operator*(const qsim::math::composition<Obj...>& comp, const V& v) {
    return comp.linear(v, qsim::helper::application);
}*/

