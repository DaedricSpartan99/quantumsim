#include "bounded_particle.hpp"
#include "debug.hpp"

using namespace qsim2d;


const ScalarField BoundedParticle::natural_potential(const ScalarField& V, double hbar) {

  return [V, hbar](vertex_t x) -> double {
    return V(x/hbar);
  };
}


// Most general constructor
BoundedParticle::BoundedParticle(std::weak_ptr<IslandMesh> mesh, 
                  const double mass, 
                  const cpx_vector& initial_state, 
                  const ScalarField& initial_potential, 
                  const double hbar,
                  std::shared_ptr<const Integrator> integrator) 
  : 
    mesh(mesh),
    state(initial_state),
    mass(mass),
    hbar(hbar),
    laplace(mesh, unit_function<double>, integrator),
    laplace_matrix(laplace.generate_matrix() * (- hbar * hbar / (2*mass))),
    unit_mass(mesh, unit_function<double>, integrator),
    unit_mass_matrix(unit_mass.generate_matrix()),
    potential_mass(mesh, initial_potential, integrator),
    potential_mass_matrix(potential_mass.generate_matrix()),
    position_x_mass(mesh, [](vertex_t v) -> double { return v[0]; }, integrator),
    position_y_mass(mesh, [](vertex_t v) -> double { return v[1]; }, integrator),
    position_x_matrix(position_x_mass.generate_matrix()),
    position_y_matrix(position_y_mass.generate_matrix())
{
  // TODO: check initial state is normalized (squared norm)

  // TODO: check coherence of potential with the state (eig. problem)
}

// evolve changing potential
void BoundedParticle::_manage_evolution(double dt, const matrix& new_potential_mass_matrix)
{
  using namespace std::complex_literals;

  // complex half step
  complex half_step = dt / (2 * hbar * 1i);

  // create backwards euler evolution operator (new potential)
  cpx_matrix U_backwards = unit_mass_matrix - (new_potential_mass_matrix + laplace_matrix) * half_step;

  // create forward euler evolution operator (old potential)
  cpx_matrix U_forward = unit_mass_matrix + (potential_mass_matrix + laplace_matrix) * half_step;

  // solve system of equations U[-1/2] * psi[n+1] = U[1/2] * psi[n]
  cpx_vector half_state = U_forward * state;
  state = U_backwards.colPivHouseholderQr().solve(half_state);

  //npdebug("Quality of the solution: ", (U_backwards * state - half_state).norm())
}


// evolve changing potential
void BoundedParticle::advance(double dt, ScalarField potential)
{
  // change potential to the new one
  potential_mass.set_field(potential);

  // generate next matrix
  matrix potential_mass_next = potential_mass.generate_matrix();
  
  // evolve 
  _manage_evolution(dt, potential_mass_next);

  // set potential matrix
  potential_mass_matrix = std::move(potential_mass_next);
}

// evolve without changing potential
void BoundedParticle::advance(double dt)
{
  // evolve
  _manage_evolution(dt, potential_mass_matrix);
}

// TODO: set probability 1 to triangle
void collapse(vertex_t position) {
}

// access state
void BoundedParticle::set_state(const cpx_vector& state) 
{
  this->state = state;
}

const cpx_vector& BoundedParticle::get_state() const 
{
  return state;
}

// access properties
double BoundedParticle::get_energy() const 
{
  // construct hamiltonian from current potential
  matrix H = potential_mass_matrix + laplace_matrix;
  
  // return energy estimation
  complex res = state.dot(H * state);

  //npdebug("Energy: ", res)
  //npdebug("Potential norm: ", potential_mass_matrix.norm())

  return res.real();
}

vertex_t BoundedParticle::get_position() const 
{
  complex x = state.dot(position_x_matrix * state);
  complex y = state.dot(position_y_matrix * state);

  //npdebug("Position X: ", x)
  //npdebug("Position Y: ", y)

  return vertex_t{x.real(), y.real()};
}

// get probability density: probability for each vertex
vector BoundedParticle::get_prob_density() const 
{
}

double BoundedParticle::get_norm() const {
  
  complex res = state.dot(unit_mass_matrix * state);

  return res.real();
}


