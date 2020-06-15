#include "grid/qsystem2D.hpp"
#include "constants.hpp"
#include "potential.hpp"

#include "grid/integrators.hpp"

#include "debug.hpp"

using namespace qsim::grid;
using namespace qsim::math;

const diagonals<wave_t, 3> qsystem2D::A = diagonals<wave_t, 3>({sdiag_entry(-1, 1.0), sdiag_entry(0, -2.0), sdiag_entry(1, 1.0)});
const diagonals<wave_t, 2> qsystem2D::A_P = diagonals<wave_t, 2>({sdiag_entry(-1, -1.0), sdiag_entry(1, 1.0)});

// init pack
wave_grid qsystem2D::init_pack::generate(double dx, double dy) const {
    matrix<wave_t> psi(N, M, wave_t(0.0));

    // construct it using the analytic expression
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j)
            psi.at(i,j) = f(i * dx, j * dy);
    }

    return psi;
}

// constructor
qsystem2D::qsystem2D(double _m, 
                     double _dx, double _dy,
                     std::shared_ptr<potential<size_t,size_t>> _V,
                     const init_pack& init,
                     std::shared_ptr<evolver> _evolver,
                     double hbar
                    )
    : qbi_gridsystem(_m, init.generate(_dx,_dy), _V, _evolver, hbar),
      dx(_dx), dy(_dy)
{
    normalize();
}

// allow to set boundaries
void qsystem2D::evolve(double dt) {
    wave = m_evolver->evolve(*this, dt);
}

diagonals<wave_t, 3> qsystem2D::H_zero_x() const {
    return wave_t(-pow(hbar() / dx, 2) / (2 * mass())) * A;
}

diagonals<wave_t, 3> qsystem2D::H_zero_y() const {
    return wave_t(-pow(hbar() / dy, 2) / (2 * mass())) * A;
}

diagonals<wave_t, 2> qsystem2D::Px() const {
    using namespace std::complex_literals;
    return (-1i * hbar() / (2.0 * dx)) * A_P;
}

diagonals<wave_t, 2> qsystem2D::Py() const {
    using namespace std::complex_literals;
    return (-1i * hbar() / (2.0 * dy)) * A_P;
}

diag_functor<wave_t> qsystem2D::potential_on_row(size_t i) const {
    return std::function<wave_t (size_t)>([&] (size_t k) -> double { return V()(i,k); });
}

diag_functor<wave_t> qsystem2D::potential_on_column(size_t j) const {
    return std::function<wave_t (size_t)>([&] (size_t k) -> double { return V()(k,j); });
}


// normalize the wave function
double qsystem2D::norm() const {
    double out(0);
    for (size_t i = 0; i < N(); ++i) {
        for (size_t j = 0; j < M(); ++j) {
            out += std::norm(wave(i,j));
            //npdebug("norm accumulation: i = ", i, ", j = ", j, ", value = ", out)
        }
    }

    npdebug("norm = ", out * dx * dy)

    return out * dx * dy;
}

std::pair<double,double> qsystem2D::position() const {
    std::pair<double,double> out(0, 0);
    for (size_t i = 0; i < N(); ++i) {
        for (size_t j = 0; j < M(); ++j) {
            double val = std::norm(wave(i,j));
            out.first += x(i) * val;
            out.second += y(j) * val;
        }
    }

    out.first *= dx * dy;
    out.second *= dx * dy;
    return out;
}

std::pair<double,double> qsystem2D::momentum() const {
    qsim::wave_t px(0), py(0);
    wave_grid& wav = const_cast<wave_grid&>(wave);

    for (size_t i = 0; i < N(); ++i)
        px += std::conj(wav.get_row(i)) * (Px() * wav.get_row(i));

    for (size_t j = 0; j < M(); ++j)
        py += std::conj(wav.get_column(j)) * (Py() * wav.get_column(j));

    using namespace std::complex_literals;

    //px *= - 1i * hbar() * dy / 2.0;
    //py *= - 1i * hbar() * dx / 2.0;
    px *= dx * dy;
    py *= dx * dy;

    if (abs(px.imag()) > qsim::machine_prec || abs(py.imag()) > qsim::machine_prec) {
        //npdebug("Values: (", px, ", ", py, ")")
        throw std::runtime_error("Momentum computation isn't fully real");
    }

    return std::pair<double, double>(px.real(), py.real());
}


// implementations
double qsystem2D::energy() const {
    qsim::wave_t E(0);
    wave_grid& wav = const_cast<wave_grid&>(wave);

    //npdebug("Constant wave: ", &wave)
    //npdebug("Non-constant wave: ", &wav)
    for (size_t i = 0; i < N(); ++i) {
        auto conj = std::conj(wav.get_row(i));
        E += conj * (H_zero_x() * wav.get_row(i));
        E += conj * (potential_on_row(i) * wav.get_row(i));
    }

    for (size_t j = 0; j < M(); ++j) {
        auto conj = std::conj(wav.get_column(j));
        E += conj * (H_zero_y() * wav.get_column(j));
        npdebug("Operating on column (phase 3) j = ", j, ", A = ", (wav.get_column(j) * wav.get_column(j)))
        E += conj * (potential_on_column(j) * wav.get_column(j));
        npdebug("Operating on column (phase 4) j = ", j, ", A = ", (wav.get_column(j) * wav.get_column(j)))
    }

    if (abs(E.imag()) > qsim::machine_prec) {
        //npdebug("Value: ", E)
        throw std::runtime_error("Energy computation isn't fully real");
    }

    return E.real();
}

void qsystem2D::set_delta_x(double _dx) {
    if (_dx <= 0)
        throw std::invalid_argument("dx must be positive");
    dx = _dx;
}

double qsystem2D::delta_x() const {
    return dx;
}

void qsystem2D::set_delta_y(double _dy) {
    if (_dy <= 0)
        throw std::invalid_argument("dy must be positive");
    dy = _dy;
}

double qsystem2D::delta_y() const {
    return dy;
}
 
