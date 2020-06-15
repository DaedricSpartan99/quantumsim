#include "grid/qsystem1D.hpp"
#include "constants.hpp"
#include "potential.hpp"

#include "debug.hpp"

using namespace qsim::grid;
using namespace qsim::math;

const diagonals<wave_t, 3> qsystem1D::A = diagonals<wave_t, 3>({sdiag_entry(-1, 1.0), sdiag_entry(0, -2.0), sdiag_entry(1, 1.0)});
const diagonals<wave_t, 2> qsystem1D::P_templ = diagonals<wave_t, 2>({sdiag_entry(-1, -1.0), sdiag_entry(1, 1.0)});

qsystem1D::qsystem1D(double _m, 
                     double _dx,
                     std::shared_ptr<potential<size_t>> _V,
                     const init_pack& init,
                     std::shared_ptr<evolver> _evolver,
                     double _hbar
                     ) : qgridsystem(_m, init.generate(_dx), _V, _evolver, _hbar),
                        dx(_dx)
                        // compose the hemiltonian matrix
    {
        normalize();
    }


wave_vector qsystem1D::init_pack::generate(double dx) const {
    wave_vector w(N, qsim::wave_t(0.0));

    // construct it using the analytic expression
    for (size_t i = 0; i < N; ++i) {
        w[i] = f(i * dx);
    }

    return w;
}

diagonals<wave_t, 3> qsystem1D::H_zero() const {
    wave_t h0 = - pow(hbar()/dx, 2) / (2 * mass());
    diagonals<wave_t,3> out(h0 * A);
    return out;
}

H_matrix_1D qsystem1D::H() const {
    // generated in-place because it's light-weight
    return H_matrix_1D(
                wave_t(0.0), 
                H_zero(), 
                std::function<wave_t (size_t)>([&] (size_t k) -> double { return V()(k); }) 
                      );
}

P_matrix_1D qsystem1D::P() const {
    using namespace std::complex_literals;
    return (-1i * hbar()) * P_templ;
}

void qsystem1D::evolve(double dt) {
    wave = m_evolver->evolve(*this, dt);
}

double qsystem1D::norm() const {
    return wave.square_norm() * dx;
}

void qsystem1D::replace_wave(const init_pack& init) {
    wave = init.generate(dx);
    normalize();
}

void qsystem1D::set_delta(double _dx) {
    if (_dx <= 0)
        throw std::invalid_argument("dx must be positive");
    dx = _dx;
}

double qsystem1D::delta() const {
    return dx;
}

double qsystem1D::energy() const {

    wave_t result = wave * (H() * wave); 

    if (abs(result.imag()) > qsim::machine_prec) {
        npdebug("psi: ", wave)
        npdebug("H * psi: ", (H() * wave))
        npdebug("Complex part of the energy: ", result.imag())
        throw result; // TODO, a real error
    }

    return result.real() * dx;
}

double qsystem1D::position() const {
    double pos(0);

    for (size_t i = 0; i < wave.size(); ++i)
        pos += x(i) * std::norm(wave[i]);

    return pos *= dx;
}

double qsystem1D::momentum() const {
    wave_t mom = wave * (P() * wave);

    if (abs(mom.imag()) > qsim::machine_prec)
        throw mom; 

    return mom.real() * dx;
}

double qsystem1D::probability(size_t beg, size_t end) const {
    if (end < beg)
        std::swap(beg, end);

    if (beg >= wave.size())
        return 0.0;

    if (end > wave.size())
        end = wave.size();

    double pos(0);

    for (size_t i = beg; i < end; ++i)
        pos += std::norm(wave[i]);

    return pos *= dx;
}

wave_vector::iterator qsystem1D::begin() {
    return wave.begin();
}

wave_vector::iterator qsystem1D::end() {
    return wave.end();
}

wave_vector::const_iterator qsystem1D::begin() const {
    return wave.begin();
}

wave_vector::const_iterator qsystem1D::end() const {
    return wave.end();
}
