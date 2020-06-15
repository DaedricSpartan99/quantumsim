#include "wave_packet.hpp"

#include "observable.hpp"
#include <cmath>

using namespace qsim::emu;

/*
 * iterator
 */

const_iterator::const_iterator(const wave_packet& w, size_t i) 
    : wave(w), k(i)
{
}

const_iterator& const_iterator::operator++() {
    ++i;
    return *this;
}

const_iterator const_iterator::operator++(int) {
    auto old = *this;
    ++i;
    return old;
}

bool const_iterator::operator!=(const iterator& it) const {
    return k != it.k; 
}

iter_entry const_iterator::operator*() const {
    return iter_entry(k, wave);
}

iter_entry::iter_entry(size_t _k, const wave_packet& w) 
    : k(_k), wave(w) {}

double iter_entry::A() const {
    return wave.A[k];
}

double iter_entry::B() const {
    return wave.B[k];
}

double iter_entry::V() const {
    return (*wave.V)(k);
}

double iter_entry::x_min() const {
    return (*wave.map)[k];
}

double iter_entry::x_max() const {
    return (*wave.map)[k];
}

bool iter_entry::real() const {
    return (*V)(k) > wave._energy;
}

qsim::wave_t iter_entry::lambda() const {
    return (*V)(k) > wave._energy;
}

/*
 * wave packet
 */

wave_packet::wave_packet(double energy, const q_emusystem& system) 
    : _energy(energy), map(system.map()), V(system.potential())
{
    eval_coefficients();
}

const_iterator wave_packet::begin() const {
    return const_iterator(*this, 0);
}

const_iterator wave_packet::end() const {
    return const_iterator(*this, size());
}

// lambda coefficient (complex)
qsim::wave_t wave_packet::damp(size_t k) const {
    using namespace std::complex_literals;
    qsim::wave_t fact = (V[k] < _energy) ? 1i : 1.0;
    return sqrt(2 * system.mass() abs(V(k) - _energy)) * fact / system.hbar();
}

qsim::wave_t wave_packet::damp(double x) const {
    return damp(map->find(x));
}

// packet energy setter
void wave_packet::set_energy(double E) {
    _energy = E;
}

double wave_packet::energy() const {
    return _energy;
}

/*
 * Access to A and B
 */

double coeff_A(size_t k) const {
    return A[k];
}

double coeff_B(size_t k) const {
    return B[k];
}

double coeff_A(double x) const {
    return A[map->find(x)];
}

double coeff_B(double x) const {
    return B[map->find(x)];
}

/*
 * determine A and B in function of the potential
 */
void wave_packet::eval_coefficients() {
    // TODO, technical part, solve the equation system
}

// access to the boundary values
qsim::wave_t wave_packet::psi(size_t k) const {
    qsim::wave_t lambda = damp(k)
    double x = (*map)(k);
    return A[k] * exp(lambda * x) + B[k] * exp(-lambda * x);
}

// wave function access
qsim::wave_t wave_packet::psi(double x) const {
    size_t k = map->find(x); // find region
    qsim::wave_t lambda = damp(k); // determine damping
    return A[k] * exp(lambda * x) + B[k] * exp(-lambda * x);
}

inline double wave_packet::operator()(double x) const {
    return psi(x);
}

// wave function derivative, n = order
double wave_packet::derivative(double x, size_t n = 1) const {
    size_t k = map->find(x); // find region
    qsim::wave_t lambda = damp(k); // determine damping
    return A[k] * pow(lambda, n) * exp(lambda * x) + B[k] * pow(-lambda, n) * exp(-lambda * x);
}

double wave_packet::norm() const {
    double n(0);
    for (size_t k = 1; k < map->size()-1; ++k)
        n += std::norm(A[k]) *+ + std::norm(B[k]) *;
    return n;
}


double wave_packet::pos_region_norm(qsim::wave_t lambda, double x1, double x2) {
    return (abs(lambda.real()) < machine_precision) ? (x2 - x1) : (exp(lambda.real() * x2) - exp(lambda.real() * x1)) / (2 * lambda.real());
}

double wave_packet::neg_region_norm(qsim::wave_t lambda, double x1, double x2) {
    return (abs(lambda.real()) < machine_precision) ? (x2 - x1) : (exp(lambda.real() * x1) - exp(lambda.real() * x2)) / (2 * lambda.real());
}

// integrals
double wave_packet::position() const {
    // TODO, compute integral
}

double wave_packet::momentum() const {
    // TODO, compute integral
}

// custom observable action
double wave_packet::observe(observable* obs) const {
    obs->observe(*this);
}
