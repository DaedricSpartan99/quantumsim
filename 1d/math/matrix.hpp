#pragma once

#include <vector>
#include <limits>
#include <deque>
#include <functional>

#include "debug.hpp"

namespace qsim::math {

    template <typename T>
    class matrix;

    template <typename T>
    class square_matrix;

    /*
     * General matrix class
     */
    template <typename T>
    class basic_matrix {
    public:
        
        basic_matrix() = default;
        basic_matrix(basic_matrix&&) = default;
        basic_matrix(const basic_matrix&) = default;

        virtual ~basic_matrix() = default;

        // without boundary check
        virtual T& operator()(size_t i, size_t j) = 0;
        virtual const T& operator()(size_t i, size_t j) const = 0;

        virtual size_t rows_nb() const = 0;
        virtual size_t cols_nb() const = 0;
        
        // with boundary check
        const T& at(size_t i, size_t j) const;
        T& at(size_t i, size_t j);

        /*
         * Swap content values
         */
        void swap(basic_matrix& other);        

        /*
         * Assignment operator
         */
        basic_matrix& operator=(const basic_matrix&);
        basic_matrix& operator=(basic_matrix&&) = default;

        /*
         * Arithmetic operator overloading
         */

        basic_matrix& operator+=(const basic_matrix&);
        basic_matrix& operator-=(const basic_matrix&);

        basic_matrix& operator*=(const T&);
        basic_matrix& operator/=(const T&);

        /*
         * Logical operator overloading
         */

        bool operator==(const basic_matrix&) const;
        bool operator!=(const basic_matrix&) const;
    };
// operator forwarding macro
#define FORWARD_BASIC_MATRIX(derived, T)                       \
                                                                \
    inline derived& operator+=(const basic_matrix<T>& other) { \
        basic_matrix<T>::operator+=(other);                 \
        return *this;                                       \
    }                                                       \
                                                            \
    inline derived& operator-=(const basic_matrix<T>& other) { \
        basic_matrix<T>::operator-=(other);                 \
        return *this;                                       \
    }                                                       \
                                                            \
    inline derived& operator*=(const T& scalar) {           \
        basic_matrix<T>::operator*=(scalar);                \
        return *this;                                       \
    }                                                       \
    inline derived& operator/=(const T& scalar) {           \
        basic_matrix<T>::operator/=(scalar);                \
        return *this;                                       \
    }


    template <typename T>
    class vector_access {
        public:

            virtual ~vector_access() = default;
            
            virtual size_t size() const = 0;

            virtual const T& operator[](size_t j) const = 0;
            virtual T& operator[](size_t j) = 0;

            const T& at(size_t j) const;
            T& at(size_t j);

            /*
             * Scalar product
             */

            T operator*(const vector_access<T>&) const;
            T scalar(const vector_access<T>&) const;
    };

    /*
     * Element-wise sub-matrix wrapper
     */
    template <typename T>
    class submatrix : public basic_matrix<T> {
    private:
        friend class matrix<T>;

        // reference to the matrix
        matrix<T>* ref;
        bool standalone;

    protected:
        // {begin index, size}
        const std::pair<size_t, size_t> rows, cols;

        // reference, {begin, end}, {begin, end}
        // construct a non-standalone submatrix
        submatrix(matrix<T>* _ref, std::pair<size_t, size_t> _rows, std::pair<size_t, size_t> _cols);
        
    public:
        submatrix() = delete;
        
        // move a submatrix and invalidate the other
        submatrix(submatrix&&);
        
        // construct a standalone copy
        // notice, this is the only way to construct a standalone submatrix
        submatrix(const submatrix&);

        // copy the content
        submatrix& operator=(const submatrix&);

        // move this object
        submatrix& operator=(submatrix&&) = delete;

        // delete the instance if standalone
        virtual ~submatrix();

        /*
         * target index assigning
         */
        
        virtual T& operator()(size_t i, size_t j) override;
        virtual const T& operator()(size_t i, size_t j) const override;

        virtual size_t rows_nb() const override {
            return rows.second;
        }

        virtual size_t cols_nb() const override {
            return cols.second;
        }
        
        /*
         * implicit convertion to a matrix
         */
        //operator matrix<T>() const;

        // operator forwarding
        FORWARD_BASIC_MATRIX(submatrix, T)
    };

    template <typename T>
    class row_vector : public vector_access<T>, public submatrix<T> {
            friend class matrix<T>;

            row_vector(matrix<T>* _ref, size_t row, std::pair<size_t, size_t> _cols)
                : submatrix<T>(_ref, std::pair(row, 1u), _cols) {}

            using submatrix<T>::rows_nb;
            using submatrix<T>::cols_nb;
            using submatrix<T>::at;
            using submatrix<T>::operator();

            template <class R, typename S, typename>
            friend R convert(R A, const std::function<S (const S&)>& operation);

        public:

            //row_vector(const row_vector& other) : submatrix<T>(other) {}
            //row_vector(row_vector&& other) : submatrix<T>(std::forward<row_vector>(other)) {}

            // TODO hope it calls superclass
            //row_vector(const row_vector&) = default;
            //row_vector(row_vector&&) = default;

            using vector_access<T>::at;
            
            // TODO hope it calls superclass
            /*inline row_vector& operator=(const row_vector& other) {
                submatrix<T>::operator=(other);
                return *this;
            }*/

            /*row_vector& operator=(row_vector&& other) {
                submatrix<T>::operator=(other);
                return *this;
            }*/

            //operator matrix<T>() const;
            
            virtual size_t size() const override {
                return cols_nb();
            } 

            virtual const T& operator[](size_t j) const override {
                return (*this)(0, j);
            } 

            virtual T& operator[](size_t j) override {
                return (*this)(0, j);
            } 

            // operator forwarding
            FORWARD_BASIC_MATRIX(row_vector, T)
    };

    template <typename T>
    class column_vector : public vector_access<T>, public submatrix<T> {
            friend class matrix<T>;

            column_vector(matrix<T>* _ref, size_t col, std::pair<size_t, size_t> _rows)
                : submatrix<T>(_ref, _rows, std::pair(col, 1u)) {}

            using submatrix<T>::rows_nb;
            using submatrix<T>::cols_nb;
            using submatrix<T>::at;
            using submatrix<T>::operator();
            
            template <class R, typename S, typename>
            friend R convert(R A, const std::function<S (const S&)>& operation);

        public:

            //column_vector(const column_vector& other) : submatrix<T>(other) {}
            //column_vector(column_vector&& other) : submatrix<T>(std::forward<column_vector>(other)) {}

            // TODO hope it calls superclass
            //column_vector(const column_vector&) = default;
            //column_vector(column_vector&&) = default;

            /*inline column_vector& operator=(const column_vector& other) {
                submatrix<T>::operator=(other);
                return *this;
            }

            inline column_vector& operator=(column_vector&& other) {
                submatrix<T>::operator=(other);
                return *this;
            }*/

            using vector_access<T>::at;

            //operator matrix<T>() const;
            
            virtual size_t size() const override {
                return rows_nb();
            } 

            virtual const T& operator[](size_t i) const override {
                return (*this)(i, 0);
            } 

            virtual T& operator[](size_t i) override {
                return (*this)(i, 0);
            } 

            // operator forwarding
            FORWARD_BASIC_MATRIX(column_vector, T)
    };
    
    template <typename T>
    class table_row : private std::vector<T>, public basic_matrix<T>, public vector_access<T> {
        friend class matrix<T>;
        friend class square_matrix<T>;

        virtual T& operator()(size_t i, size_t) override {
            return (*this)[i];
        }

        virtual const T& operator()(size_t i, size_t) const {
            return (*this)[i];
        }

        virtual size_t rows_nb() const {
            return 1u;
        }

        virtual size_t cols_nb() const {
            return this->size();
        }

        template <class R, typename S, typename>
        friend R convert(R A, const std::function<S (const S&)>& operation);

    public:
        using std::vector<T>::vector;

        virtual T& operator[](size_t i) override {
            return std::vector<T>::operator[](i);
        }

        virtual const T& operator[](size_t i) const {
            return std::vector<T>::operator[](i);
        }

        virtual size_t size() const override {
            return std::vector<T>::size();
        }

        // operator forwarding
        FORWARD_BASIC_MATRIX(table_row, T)
    };

    /*
     * General matrix class, dynamically allocated
     */
    template <typename T>
    class matrix : public basic_matrix<T> {
    protected:

        struct table_t : public std::vector<table_row<T>> {
            using std::vector<table_row<T>>::vector;
            using std::vector<table_row<T>>::push_back;
            using std::vector<table_row<T>>::pop_back;
            using std::vector<table_row<T>>::reserve;
        };
    
        // rows table
        table_t table;
        
        // only if needed
        matrix() = default;

    public:

        matrix(size_t N, size_t M, T init = 0);

        matrix(std::initializer_list<std::initializer_list<T>> list);
        
        //matrix(const matrix&) = default;

        // move table
        //matrix(matrix&&) = default;
        //matrix& operator=(matrix&&) = default;
        //matrix& operator=(const matrix&) = default;

        // construct copying a submatrix
        matrix(const submatrix<T>& sub);

        virtual ~matrix() = default;
        
        // no, boundary check
        virtual const T& operator()(size_t i, size_t j) const override {
            return table[i][j];
        }

        virtual T& operator()(size_t i, size_t j) override {
            return table[i][j];
        }

        virtual size_t rows_nb() const override {
            return table.size();
        }

        virtual size_t cols_nb() const override {
            return table[0].size();
        }
        
        /*
         * Create a sub-matrix with the given indexes
         */
        submatrix<T> restrict(std::pair<size_t,size_t> rows = all, std::pair<size_t, size_t> cols = all); 

        /*
         * Create a sub-matrix with a column indexing restriction
         */
        row_vector<T> get_row(size_t i, std::pair<size_t, size_t> cols);
        
        /*
         * Row descriptor without a column restriction
         */
        inline table_row<T>& get_row(size_t i) {
            return table[i];
        }

        inline const table_row<T>& get_row(size_t i) const {
            return table[i];
        }

        /*
         * Create a sub-matrix with a row indexing restriction
         */
        column_vector<T> get_column(size_t j, std::pair<size_t,size_t> rows = all);

        /*
         * swap the matrix internally
         */

        void swap_rows(std::pair<size_t, size_t> swp);
        
        // using a column range
        void swap_rows(std::pair<size_t, size_t> swp, std::pair<size_t, size_t> col_range);

        /*
         * access aliases
         */
        
        inline submatrix<T> operator()(std::pair<size_t,size_t> rows, std::pair<size_t,size_t> cols) {
            return restrict(rows, cols);
        }
        
        /*
         * helpers for restrict
         */
        static constexpr std::pair<size_t, size_t> all = {0, 0};

        // return a single line 
        static constexpr std::pair<size_t, size_t> single(size_t k) {
            return std::pair<size_t,size_t>(k, k+1);
        }
        
        // operator forwarding
        FORWARD_BASIC_MATRIX(matrix, T)
    };

    /*
     * Square matrix specialization
     */
    template<typename T>
    class square_matrix : public matrix<T> {
    public:
        square_matrix(size_t N, T init = 0);

        square_matrix(std::initializer_list<std::initializer_list<T>> list);

        //square_matrix(const square_matrix&) = default;
        //square_matrix(square_matrix&&) = default;

        //square_matrix& operator=(const square_matrix&) = default;
        //square_matrix& operator=(square_matrix&&) = default;

        inline size_t size() const {
            return matrix<T>::rows_nb();
        }
        
        /*
         * Generate the identity matrix
         */
        static square_matrix eye(size_t N);

        // operator forwarding
        FORWARD_BASIC_MATRIX(square_matrix, T)
    };
    
    // convert using an element-wise operation
    template <class R, typename T, typename = typename std::enable_if<std::is_base_of<basic_matrix<T>, R>::value>::type >
    R convert(R A, const std::function<T (const T&)>& operation);
}

/*
 * std::abs overload
 */

namespace std {


    /*
     * Maximal indexes corresponding to the maximal value
     */
    template<typename T>
    std::pair<size_t, size_t> max(const qsim::math::basic_matrix<T>& mat);

    /*
     * Element-wise absolute value
     */

    // element-wise absolute value
    template <typename T>
    qsim::math::matrix<double> abs(qsim::math::matrix<T> A);

    // complex case adjoint absolute value
    template <typename T>
    qsim::math::submatrix<double> abs(qsim::math::submatrix<T> A);

    // complex case adjoint absolute value
    template <typename T>
    qsim::math::row_vector<double> abs(qsim::math::row_vector<T> A);

    // complex case adjoint absolute value
    template <typename T>
    qsim::math::column_vector<double> abs(qsim::math::column_vector<T> A);

    // complex case adjoint absolute value
    template <typename T>
    qsim::math::table_row<double> abs(qsim::math::table_row<T> A);

    /*
     * Element-wise complex conjugate
     */

    // complex case adjoint absolute value
    template <typename T>
    qsim::math::matrix<std::complex<T>> conj(qsim::math::matrix<std::complex<T>> A);

    // complex case adjoint absolute value
    template <typename T>
    qsim::math::submatrix<std::complex<T>> conj(qsim::math::submatrix<std::complex<T>> A);

    // complex case adjoint absolute value
    template <typename T>
    qsim::math::row_vector<std::complex<T>> conj(qsim::math::row_vector<std::complex<T>> A);

    // complex case adjoint absolute value
    template <typename T>
    qsim::math::column_vector<std::complex<T>> conj(qsim::math::column_vector<std::complex<T>> A);

    // complex case adjoint absolute value
    template <typename T>
    qsim::math::table_row<std::complex<T>> conj(qsim::math::table_row<std::complex<T>> A);
}

/*
 * Operator overloading for matrix
 */

template<typename T, class Vector>
Vector operator*(const qsim::math::square_matrix<T>& A, const Vector& x);

/*
 * LU decomposition section
 */
namespace qsim::math {

    namespace helper {
        typedef std::deque<std::pair<size_t, size_t>> pivot_t;
        
        template<typename T, class Matrix>
        struct LU_output {
            square_matrix<T> L;
            square_matrix<T> U;

            // pivoting swap operations stack
            pivot_t P;

            LU_output(const Matrix& A)
                : L(square_matrix<T>::eye(A.size())), U(A) {}
        };

        
        // gauss elimination of lower triangular
        template<typename T, class Vector>
        Vector solve_lower(const square_matrix<T> L, Vector y);

        // gauss elimination of upper triangular
        template<typename T, class Vector>
        Vector solve_upper(const square_matrix<T> U, Vector x);
    }
    
    /*
     * Compute LU decomposition
     * A = intering square matrix
     * Requirements: size(), static_cast to square_matrix<T>
     */
    template<typename T, class Matrix>
    helper::LU_output<T, Matrix> LU_decomposition(const Matrix& A);
       
    /*
     * Solve A * s = b system
     */ 
    template<typename T, class Matrix, class Vector>
    Vector solve(const Matrix& A, Vector b);
}


namespace qsim::math {

/*
 * Operator+
 */ 
template<typename T>
const qsim::math::matrix<T> operator+(qsim::math::matrix<T>, const qsim::math::basic_matrix<T>&);

template<typename T>
const qsim::math::submatrix<T> operator+(qsim::math::submatrix<T>, const qsim::math::basic_matrix<T>&);

template<typename T>
const qsim::math::row_vector<T> operator+(qsim::math::row_vector<T>, const qsim::math::basic_matrix<T>&);

template<typename T>
const qsim::math::table_row<T> operator+(qsim::math::table_row<T>, const qsim::math::basic_matrix<T>&);

template<typename T>
const qsim::math::column_vector<T> operator+(qsim::math::column_vector<T>, const qsim::math::basic_matrix<T>&);

/*
 * Operator-
 */
template<typename T>
const qsim::math::matrix<T> operator-(qsim::math::matrix<T>, const qsim::math::basic_matrix<T>&);

template<typename T>
const qsim::math::submatrix<T> operator-(qsim::math::submatrix<T>, const qsim::math::basic_matrix<T>&);

template<typename T>
const qsim::math::row_vector<T> operator-(qsim::math::row_vector<T>, const qsim::math::basic_matrix<T>&);

template<typename T>
const qsim::math::table_row<T> operator-(qsim::math::table_row<T>, const qsim::math::basic_matrix<T>&);

template<typename T>
const qsim::math::column_vector<T> operator-(qsim::math::column_vector<T>, const qsim::math::basic_matrix<T>&);

/*
 * Operator*
 */
template<typename T>
const qsim::math::matrix<T> operator*(qsim::math::matrix<T>, const T&);

template<typename T>
const qsim::math::matrix<T> operator*(const T&, qsim::math::matrix<T>);


template<typename T>
const qsim::math::submatrix<T> operator*(qsim::math::submatrix<T>, const T&);

template<typename T>
const qsim::math::submatrix<T> operator*(const T&, qsim::math::submatrix<T>);


template<typename T>
const qsim::math::row_vector<T> operator*(qsim::math::row_vector<T>, const T&);

template<typename T>
const qsim::math::row_vector<T> operator*(const T&, qsim::math::row_vector<T>);


template<typename T>
const qsim::math::table_row<T> operator*(qsim::math::table_row<T>, const T&);

template<typename T>
const qsim::math::table_row<T> operator*(const T&, qsim::math::table_row<T>);


template<typename T>
const qsim::math::column_vector<T> operator*(qsim::math::column_vector<T>, const T&);

template<typename T>
const qsim::math::column_vector<T> operator*(const T&, qsim::math::column_vector<T>);

/*
 * Operator/
 */
template<typename T>
const qsim::math::matrix<T> operator/(qsim::math::matrix<T>, const T&);

template<typename T>
const qsim::math::submatrix<T> operator/(qsim::math::submatrix<T>, const T&);

template<typename T>
const qsim::math::row_vector<T> operator/(qsim::math::row_vector<T>, const T&);

template<typename T>
const qsim::math::table_row<T> operator/(qsim::math::table_row<T>, const T&);

template<typename T>
const qsim::math::column_vector<T> operator/(qsim::math::column_vector<T>, const T&);

}

// implementation
#include "matrix_template.hpp"

