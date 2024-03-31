#include "load.hpp"

#include "debug.hpp"

using namespace qsim2d;

std::array<std::function<double(const vertex_t&)>, 3> Load::basis = {
        [](const vertex_t& z) -> double { return 1 - z[0] - z[1]; }, 
        [](const vertex_t& z) -> double { return z[0]; },
        [](const vertex_t& z) -> double { return z[1]; }
      };

Load::Load(
    std::weak_ptr<const AbstractMesh> mesh, 
    const ScalarField& field,
    std::shared_ptr<const Integrator> integrator
    ) 
  : VectorComponent(), FEComponent(mesh, field), integrator(integrator)
{
}

vector Load::generate_vector() const {

  auto mesh_ptr = mesh.lock();

  const index_t N_active_vertex = mesh_ptr->n_active_vertices();
  const index_t N_triangles = mesh_ptr->n_triangles();

  // initialize with a null vector
  vector l = vector::Zero(N_active_vertex);

  for (index_t k = 0; k < N_triangles; ++k) {

    // get triangle structure
    const triangle_t& K = mesh.lock()->get_triangle(k);

    // get pre-computed metric structure
    const auto& metric = metrics[k];

    // store triangle contribution
    Eigen::Vector3d l_k;

    for(index_t i = 0; i < 3; ++i) {

        // integrate field over basis i
        double field_integral = integrator->integrate([&](vertex_t z) -> double {
          return field(metric.transform(z)) * basis[i](z);
        }); 


        // update local matrix
        l_k(i) = metric.abs_detB * field_integral;
    }

    // update global matrix
    for(index_t i = 0; i < 3; ++i) {

      // determine target in the edge
      const index_t N_i = K[i];

      // if N_i is out of scope: do nothing
      if (N_i >= N_active_vertex) {
        continue;
      }

      l(N_i) += l_k(i);
    }

  } 
  
  return l;
}

