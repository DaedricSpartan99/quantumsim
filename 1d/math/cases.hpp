#pragma once

#include <initializer_list>
#include "matrix.hpp"

namespace qsim::math {

    template<typename T>
    class cases : protected std::map<size_t, T> {
    protected:
        const T def_val;
    public:
        virtual T& get(size_t k) = 0;
        virtual const T& get(size_t k) const = 0;

        cases(std::initializer_list<std::pair<size_t,T>> list = {}) : def_val(0) {
            for (const auto& pair : list)
                emplace(pair.first, pair.second);
        }
        
        // construct a null object
        cases(T&& val) : def_val(std::move(val)) {}

        using std::map<size_t, T>::emplace;
    };
    
    template<typename T>
    class match_cases : public cases<T> {
    public:
        virtual T& get(size_t k) override {
            return std::map<size_t, T>::operator[](k);
        }

        const T& get(size_t k) const {
            auto it = std::map<size_t, T>::find(k);
            return (it != std::map<size_t, T>::end()) ? *it : def_val;
        }
    };

    template<typename T>
    class interval_cases : public cases<T> {
    public:

        virtual T& get(size_t k) override {
            auto it = std::map<size_t, T>::lower_bound(k);
            if (it != std::map<size_t, T>::end())
                return *it;
            else
                throw std::out_of_range("(interval_case<T>::get) element not found");
        }

        const T& get(size_t k) const {
            auto it = std::map<size_t, T>::lower_bound(k);
            return (it != std::map<size_t, T>::end()) ? *it : def_val;
        }
    };

    template<typename T>
    class cases_square_matrix : public interval_cases<match_cases<T>> {

        const size_t N;

    public:
        cases_square_matrix(size_t _N, std::initializer_list<std::pair<cases_entry*, cases<T>>> list = {})
            : cases(list), N(_N) {}

        inline size_t size() const {
            return N;
        }

        inline T& at(size_t i, size_t j) {
            return get(i).get(j);
        }

        inline const T& at(size_t i, size_t j) const {
            return get(i).get(j);
        }
        
        /*
         * Construct a square_matrix using this definition
         */
        operator square_matrix<T>() const {
            square_matrix<T> out(N);

            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    out(i,j) = at(i, j);
                }
            }

            return out;
        }
    };
}
