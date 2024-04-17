#include "types.hpp"
#include "island_mesh.hpp"
#include <iostream>
#include <vector>
#include "solvers/bounded_particle.hpp"
#include <fstream>
#include <cmath>

#include "debug.hpp"

using namespace qsim2d;

std::shared_ptr<IslandMesh> generate_mesh(int N);
cpx_vector initialize_state(std::shared_ptr<IslandMesh>, std::function<complex(vertex_t)>);

int main() {

  // mesh parameters
  const int N = 10;
  const double h = 1. / N;

  // physical parameters
  const double mass = 1.;
  const int wave_n = 1;
  const int wave_m = 1;
  const double E = M_PI * M_PI * (wave_n * wave_n + wave_m * wave_m) / 2. ;

  std::cout << "Expected energy: " << E << std::endl;
  
  // time division parameters
  const int Nt = 50;
  const double dt = 1. / Nt;

  // define analytic solution
  std::function<complex(vertex_t, double)> analytic_wave_fct = [wave_n, wave_m, E](vertex_t v, double t) -> complex {

    using namespace std::complex_literals;

    return 2 * 
      sin(wave_n * M_PI * v[0]) * 
      sin(wave_m * M_PI * v[1]) * 
      exp(-1i * E * t);
  };

  // initial wave function 
  std::function<complex(vertex_t)> initial_wave_fct = [wave_n, wave_m, analytic_wave_fct](vertex_t v) -> complex {
    return analytic_wave_fct(v, 0);
  };

  // construct island mesh
  std::shared_ptr<IslandMesh> island_mesh = generate_mesh(N);
  
  // define potential
  auto potential = [](vertex_t v) -> double { return 0; };

  // define initial state
  cpx_vector init_state = initialize_state(island_mesh, initial_wave_fct);

  // construct simulator
  BoundedParticle simulator(island_mesh, 
                            mass,
                            init_state,
                            potential);

  // fresh simulation  
  for (int n = 0; n < Nt; ++n) { 

    // step 
    std::cout << "Step: " << n << ", t = " << (n * dt) << std::endl;
    
    // check energy and normalization
    std::cout << "Energy: " << simulator.get_energy() << std::endl;
    std::cout << "Normalization: " << simulator.get_norm() << std::endl;

    // get position
    vertex_t position = simulator.get_position();

    std::cout << "<X> : " << position[0] << std::endl;
    std::cout << "<Y> : " << position[1] << std::endl;
    
    auto current_analytic = [analytic_wave_fct, n, dt](vertex_t v) -> complex {
      return analytic_wave_fct(v, n * dt);
    };

    std::cout << "Difference from analytic solution: " << (simulator.get_state() - initialize_state(island_mesh, current_analytic)).norm() << std::endl;

    std::cout << std::endl;

    // advance
    simulator.advance(dt);
  }

  std::cout << "End of simulation" << std::endl;
}

std::shared_ptr<IslandMesh> generate_mesh(int N) {

  const double h = 1. / N;

  // initialize components
  std::vector<vertex_t> vertices;
  std::vector<triangle_t> triangles;

  // draw first line of vertices 
  for (index_t j = 0; j < N+1; ++j)
    vertices.push_back(vertex_t{0., double(j) * h});
  
  // draw other vertices and bind triangles
  for (index_t i = 1; i < N+1; ++i) {

    // draw line of vertices 
    for (index_t j = 0; j < N+1; ++j)
      vertices.push_back(vertex_t{double(i) * h, double(j) * h});
    
    // construct two triangles per square (reference corner: bottom left)
    for (index_t j = 0; j < N; ++j) {

      index_t tl = (i-1) * (N+1) + j;
      index_t tr = (i-1) * (N+1) + j + 1;
      index_t bl = i * (N+1) + j;
      index_t br = i * (N+1) + j + 1;

      // bottom-left -> bottom-right -> top-right
      triangle_t T1 {bl, br, tr};

      // bottom-left -> top-left -> top-right
      triangle_t T2 {bl, tl, tr};

      // push them to collection
      triangles.push_back(T1);
      triangles.push_back(T2);
    }
  }

  // construct island mesh
  return std::make_shared<IslandMesh>(vertices, triangles);
}

cpx_vector initialize_state(std::shared_ptr<IslandMesh> mesh, std::function<complex(vertex_t)> initial_wave_fct) {

  const Mesh& internal = mesh->get_internal_mesh();

  cpx_vector state(internal.n_vertices());

  // complete with state
  for (int i = 0; i < internal.n_vertices(); ++i) {

    const vertex_t& x = internal.get_vertex(i);
    
    state[i] = initial_wave_fct(x);
  }

  return state;
}

