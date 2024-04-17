#include "mass.hpp"

using namespace qsim2d;

std::array<std::function<double(const vertex_t&)>, 3> Mass::basis = {
        [](const vertex_t& z) -> double { return 1 - z[0] - z[1]; }, 
        [](const vertex_t& z) -> double { return z[0]; },
        [](const vertex_t& z) -> double { return z[1]; }
      };



Mass::Mass(
    std::weak_ptr<const AbstractMesh> mesh, 
    const ScalarField& field,
    std::shared_ptr<const Integrator> integrator
    ) 
  : MatrixComponent(), FEComponent(mesh, field), integrator(integrator)
{
}

matrix Mass::generate_matrix() const {

  auto mesh_ptr = mesh.lock();

  const index_t N_active_vertex = mesh_ptr->n_active_vertices();
  const index_t N_triangles = mesh_ptr->n_triangles();

  // initialize with a null matrix
  matrix M = matrix::Zero(N_active_vertex, N_active_vertex);

  for (index_t k = 0; k < N_triangles; ++k) {

    // get triangle structure
    const triangle_t& K = mesh.lock()->get_triangle(k);

    // get pre-computed metric structure
    const auto& metric = metrics[k];

    // store triangle contribution
    Eigen::Matrix3d M_k;

    for(index_t i = 0; i < 3; ++i) {

      for(index_t j = 0; j < 3; ++j) {

        // integrate field * basis_i over unit triangle
        double field_integral = integrator->integrate([&](vertex_t z) -> double {
          return field(metric.transform(z)) * basis[i](z) * basis[j](z);
        }); 


        // update local matrix
        M_k(i,j) = metric.abs_detB * field_integral;
      }
    }

    // update global matrix
    for(index_t i = 0; i < 3; ++i) {

      // determine target in the edge
      const index_t N_i = K[i];

      // if N_i is out of scope: do nothing
      if (N_i >= N_active_vertex) {
        continue;
      }

      for (index_t j = 0; j < 3; ++j) {

        // determine target in the edge
        const index_t N_j = K[j];

        // if N_j is out of scope: do nothing
        if (N_j >= N_active_vertex) {
          //boundary_index.push_back(N_j);
          //boundary = true;
          continue;
        } 
        
        // assign contribution
        M(N_i, N_j) += M_k(i,j);
      }
    }

  } 
  
  return M;
}


NeumannMass::NeumannMass(
          std::weak_ptr<const Mesh> mesh, 
          const ScalarField& field,
          std::shared_ptr<const Integrator> integrator) 
  : Mass(std::static_pointer_cast<const AbstractMesh>(mesh.lock()), field, integrator)
{
}


DirichletMass::DirichletMass(
          std::weak_ptr<const IslandMesh> mesh, 
          const ScalarField& field,
          std::shared_ptr<const Integrator> integrator) 
  : Mass(std::static_pointer_cast<const AbstractMesh>(mesh.lock()), field, integrator)
{
}

