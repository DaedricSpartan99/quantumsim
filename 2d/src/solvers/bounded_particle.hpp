#pragma once

#include "abstract_solver.hpp"
#include <memory>
#include "island_mesh.hpp"
#include "stiffness.hpp"
#include "mass.hpp"

namespace qsim2d {

  /*
   * This class solves the schrodinger equation for a bounded one-particle case
   *
   * - (1/2m) laplace[u] + V * u = i * du/dt
   *
   * with dirichlet zero boundary conditions:
   *
   *  u(x) = 0 : x on boundary
   */

  class BoundedParticle {
    public:
      
      // Most general constructor
      BoundedParticle(std::weak_ptr<IslandMesh> mesh, 
                        const double mass, 
                        const cpx_vector& initial_state, 
                        const ScalarField& initial_potential,
                        const double hbar = 1.0,
                        std::shared_ptr<const Integrator> integrator = std::make_shared<const Interpolator>(GAUSS_TRI_POINTS));
    
      // evolve changing potential
      void advance(double dt, ScalarField potential);

      // evolve without changing potential
      void advance(double dt);

      // collapse after measurement
      void collapse(vertex_t position);

      // trivially access state
      void set_state(const cpx_vector&);
      const cpx_vector& get_state() const;

      // access properties
      double get_energy() const;
      vertex_t get_position() const; 

      // get probability density for each triangle in the mesh
      vector get_prob_density() const;

      // get norm: just for verification, must be 1
      double get_norm() const;

    private:

      // store mesh 
      std::weak_ptr<IslandMesh> mesh;

      // store current state
      cpx_vector state;
      
      // mass
      const double mass;

      // matrix
      const double hbar;

      // store laplace generator and matrix
      const DirichletStiffness laplace;
      const matrix laplace_matrix;

      // store mass unit M[1]
      const DirichletMass unit_mass;
      const matrix unit_mass_matrix;

      // store mass matrix of the old potential M[V_n]
      DirichletMass potential_mass;

      // store current potential mass matrix
      matrix potential_mass_matrix;

      // store position operator matrix
      const DirichletMass position_x_mass;
      const DirichletMass position_y_mass;
      const matrix position_x_matrix;
      const matrix position_y_matrix;

      /*
       * Evolver
       */

      void _manage_evolution(double dt, const matrix&);

      /*
       * Unit converters
       */

      static const ScalarField natural_potential(const ScalarField&, double hbar);
  };
}
