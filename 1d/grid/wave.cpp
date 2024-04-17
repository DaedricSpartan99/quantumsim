#include "grid/wave.hpp"
#include <algorithm>
#include <ostream>

#include "debug.hpp"

using namespace qsim::grid;

/*wave_vector::wave_vector(size_t N, wave_t value) : std::vector<wave_t>(N, value) {
    npdebug("Constructing: ", this)
}*/

wave_vector::wave_vector(const wave_vector& other) : std::vector<wave_t>(other) {
    //npdebug("Copied ", &other, " into ", this)
}

wave_vector::wave_vector(wave_vector&& other) : std::vector<wave_t>(std::forward<std::vector<wave_t>>(other)) {
    //npdebug("Moved ", &other, " into ", this, ", other empty: ", other.empty())
}

wave_vector& wave_vector::operator=(const wave_vector& other) {
    std::vector<wave_t>::operator=(other);
    return *this;
}

wave_vector::~wave_vector() {
    //npdebug("Destroy: ", this)
}

wave_vector& wave_vector::operator+=(const wave_vector& other) {
    
    const size_t N = std::min(size(), other.size()); 

    for (size_t k = 0; k < N; ++k)
        (*this)[k] += other[k];
    
    // push new values
    for (size_t k = N; other.size() > size(); ++k)
        push_back(other[k]);

    return *this;
}

wave_vector& wave_vector::operator-=(const wave_vector& other) {
    const size_t N = std::min(size(), other.size()); 

    for (size_t k = 0; k < N; ++k)
        (*this)[k] -= other[k];
    
    // push new values
    for (size_t k = N; other.size() > size(); ++k)
        push_back(-other[k]);

    return *this;
}

wave_vector& wave_vector::operator*=(const wave_t& z) {
    for (size_t k = 0; k < size(); ++k)
        (*this)[k] *= z;

    return *this;
}

wave_vector& wave_vector::operator/=(const wave_t& z) {
    for (size_t k = 0; k < size(); ++k)
        (*this)[k] /= z;

    return *this;
}

wave_t wave_vector::at(int i) const {
    return (i < 0 || static_cast<size_t>(i) >= size()) ? 
                0 : std::vector<wave_t>::operator[](static_cast<size_t>(i));
}

// namespace grid
wave_vector qsim::grid::operator+(wave_vector w, const wave_vector& other) {
    return w += other;
}

wave_vector qsim::grid::operator-(wave_vector w, const wave_vector& other) {
    return w -= other;
}

wave_vector qsim::grid::operator*(wave_vector w, const wave_t& z) {
    return w *= z;
}

wave_vector qsim::grid::operator/(wave_vector w, const wave_t& z) {
    return w /= z;
}

wave_vector qsim::grid::operator*(const wave_t& z, wave_vector w) {
    return w *= z;
}

wave_vector qsim::grid::operator/(const wave_t& z, wave_vector w) {
    return w /= z;
}


std::ostream& operator<<(std::ostream& os, const wave_vector& v) {
    using namespace std;

    os << "(";
    for (size_t k = 0; k < v.size() - 1; ++k)
        os << v[k] << ", ";
    
    if (v.size() > 0)
        os << v.back();

    os << ")";
    return os;
}

qsim::grid::wave_vector std::conj(qsim::grid::wave_vector v) {
    for (auto& val : v)
        val = std::conj(val);

    return v;
}

qsim::grid::wave_t qsim::grid::operator*(qsim::grid::wave_vector x, const qsim::grid::wave_vector& y) {
    
    wave_t result(0);

    for (size_t k = 0; k < std::min(x.size(), y.size()); ++k)
        result += std::conj(x[k]) * y[k];

    return result;
}

double qsim::grid::wave_vector::square_norm() const {

    double a(0);

    for (size_t i = 0; i < size(); ++i)
        a += std::norm((*this)[i]);

    return a;
}

