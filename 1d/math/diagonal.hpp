#pragma once

#include <array>
#include <functional>

namespace qsim::math {
    
    // constant random access
    template<typename T>
    class diagonal {
    public:

        virtual ~diagonal() = default;

        virtual T operator[](size_t m) const = 0;

        T at(size_t i, size_t j) const {
            return (i == j) ? operator[](i) : T(0);
        }
    };

    template<typename T, size_t N>
    class diag_array : public diagonal<T> {

        // array
        std::array<T,N> data;
    
    public:
        diag_array(const std::array<T,N>& _data) : data(_data) {}

        virtual T operator[](size_t m) const override {
            return data[m];
        }

        diag_array& operator*=(const T& g) {

            for (auto& elem : data)
                elem *= g;

            return *this;
        }
    };

    template<typename T>
    class diag_functor : public diagonal<T> {
        
        // functor T f(size_t)
        std::function<T (size_t)> f; 
        T gain;
    
    public:
        diag_functor(const std::function<T (size_t)>& _f, T g = 1) : f(_f), gain(g) {}

        virtual T operator[](size_t m) const override {
            return gain * f(m);
        }

        diag_functor& operator*=(const T& g) {

            gain *= g;
            return *this;
        }
    };
}

/*
 * Scalar multiplication
 */
namespace qsim::math {

template<typename T>
const qsim::math::diag_functor<T> operator*(qsim::math::diag_functor<T> A, const T& g) {
    return A *= g;
}

template<typename T>
const qsim::math::diag_functor<T> operator*(const T& g, qsim::math::diag_functor<T> A) {
    return A *= g;
}

template<typename T, size_t N>
const qsim::math::diag_array<T,N> operator*(qsim::math::diag_array<T,N> A, const T& g) {
    return A *= g;
}

template<typename T, size_t N>
const qsim::math::diag_array<T,N> operator*(const T& g, qsim::math::diag_array<T,N> A) {
    return A *= g;
}

/*
 * Matrix multiplication by diagonal
 */
template<typename T, class V>
const V operator*(const qsim::math::diagonal<T>& mat, V v) {

    for (size_t k = 0; k < v.size(); ++k)
        v[k] *= mat[k];

    return v;
}

}
