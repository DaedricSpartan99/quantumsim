#pragma once

#include <algorithm>
#include <utility>

// TODO, debug
//#include <iostream>

namespace qsim::math {

    // admit std::sort
    template <typename T> 
    struct sdiag_entry : public std::pair<long int,T> {
        
        using std::pair<long int, T>::pair;

        constexpr sdiag_entry(long int n, T x) 
            : std::pair<long int, T>(n, x) {}

        bool operator<(const sdiag_entry& other) const {
            // compare indices
            return this->first < other.first;    
        }

        bool operator<(long int index) const {
            // compare indices
            return this->first < index;    
        }

        inline long int column() const {
            return this->first;
        }

        inline const T& value() const {
            return this->second;
        }
    };
  
    /*
     * This class describes a square matrix composed by multiple
     * sub-diagonals (shifted diagonals), where
     * each diagonal admits only ONE value.
     * T = matrix field type
     * D = quantity of non-zero values in the first row
     */ 
    template <typename T, size_t D> 
    class diagonals {
        
        // avoid empty initializations
        static_assert(D > 0);

        using entry = sdiag_entry<T>;
        
        // long int = indices, 0 = diagonal
        // T = constant values at the first row
        std::array<entry,D> data;
        
        // look for a specific value via binary search O(log(D))
        auto find(long int index) const {
            auto it = std::lower_bound(data.begin(), data.end(), index);
            return (it != data.end() && it->column() == index) ? it : data.end();
        }

    public:

        diagonals(const std::array<entry,D>& values) 
            : data(values) {
                // sort by ascending index
                std::sort(data.begin(), data.end());
        }

        diagonals(const diagonals&) = default;
        diagonals(diagonals&&) = default;

        diagonals& operator=(const diagonals&) = default;
        
        // random access
        T at(size_t i, size_t j) const {
            long int k = static_cast<long int>(j) - static_cast<long int>(i);
            auto it = find(k);
            return (it != data.end()) ? it->value() : T(0);
        }

        /*
         * Output type
         */

        struct output_t : public std::pair<size_t,T> {
        
            using std::pair<size_t,T>::pair;

            inline size_t column() const {
                return this->first;
            }

            inline const T& value() const {
                return this->second;
            }
        };

        /*
         * Iterator to a row
         */
    
        struct const_iterator {

            using array_iterator = typename std::array<entry,D>::const_iterator;

            size_t row;
            array_iterator it;

            const_iterator(size_t r, array_iterator _it)
                : row(r), it(_it) {
                // skip invalid addresses
            }

            output_t operator*() {
                return output_t(row + static_cast<size_t>(it->column()), it->value());
            }

            const_iterator operator++(int) {
                auto out = *this;
                ++it;
                return out;
            }

            const_iterator operator++() {
                ++it;
                return *this;
            }

            bool operator!=(const const_iterator& other) const {
                return it != other.it;
            }
        };

        /*
         * Row descriptor, best way to access a diagonals
         * It restricts the domain to the one which is non-zero
         * N = entire size of the matrix
         */

        struct row_array {

            const diagonals<T,D>& reference;
            const const_iterator m_begin, m_end; // pre-compute the boundaries

            row_array(const diagonals<T,D>& ref, size_t row_number, size_t N) : 
                reference(ref), 
                m_begin(reference.begin(row_number)),
                m_end(reference.end(row_number, N))
            {
                // TODO, debug
                //using namespace std;
                //cout << "Begin shift: " << m_begin - reference.data.begin() << endl;
                //cout << "End shift: " << reference.data.end() - m_end << endl;
            }

            const const_iterator begin() const {
                return m_begin;
            }

            const const_iterator end() const {
                return m_end;
            }
        };

        /*
         * Row accessor and operator []
         * Advantage: a row_array is for auto loopable
         * N = entire size of the matrix
         */

        inline row_array get_row(size_t row, size_t N) const {
            return row_array(*this, row, N);
        }
        
        // multiplication by scalar O(D)
        diagonals& operator*=(const T& a) {
            for (size_t k = 0; k < data.size(); ++k)
                data[k].second *= a;

            return *this;
        }
        
        // add a multiple of the identity matrix O(log(D))
        diagonals& operator+=(const T& a) {
            auto it = this->find(0);
            if (it != data.end()) // not found
                it->second += a;
            else {
                // TODO, throw error
            }
            return *this;
        }

    private:
        /*
         * Begin, end
         */

        const_iterator begin(size_t row) const {
            auto it = data.begin();

            while(it != data.end() && (*it).column() + static_cast<long int>(row) < 0)
                ++it;

            return const_iterator(row, it);
        }

        const_iterator end(size_t row, size_t N) const {
            auto it = data.end(); // initialize to last available element

            while(it != data.begin() && (*(it-1)).column() + static_cast<long int>(row) >= static_cast<long int>(N))
                --it;

            return const_iterator(row, it);
        }
    };
    
    template<typename T, class ... Args>
    diagonals<T, sizeof...(Args) + 1> make_diagonals(sdiag_entry<T>&& first, Args&& ... args) {
        return diagonals<T, sizeof...(Args) + 1>({first, (sdiag_entry<T>(args), ...) });
    }
}

/*
 * By scalar multiplication
 */
namespace qsim::math {

template <typename T, size_t D>
const qsim::math::diagonals<T,D> operator*(qsim::math::diagonals<T,D> A, const T& a) {
    return A *= a;
}

template <typename T, size_t D> 
const qsim::math::diagonals<T,D> operator*(const T& a, qsim::math::diagonals<T,D> A) {
    return A *= a;
}

/*
 * Matrix multiplication/application
 * Enable it for all non-arithmetic types
 */
template <typename T, size_t D, class V>
const V operator*(const qsim::math::diagonals<T,D>& mat, const V& v) {
    
    V out(v); // eventually copy the size, in case it's vector like

    for (size_t m = 0; m < v.size(); ++m) {
        out[m] = 0;
        // loop through the column
        for (auto&& col : mat.get_row(m, v.size()))
            out[m] += col.value() * v[col.column()];
    }
    
    return out;
}

}
