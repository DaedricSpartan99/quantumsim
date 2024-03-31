#include "stiffness.hpp"

#include "debug.hpp"

using namespace qsim2d;

Stiffness::Stiffness(
          std::weak_ptr<const AbstractMesh> mesh, 
          const ScalarField& field,
          std::shared_ptr<const Integrator> integrator
          )
        : MatrixComponent(), FEComponent(mesh, field), integrator(integrator)
      {
      }

const Eigen::Matrix<double, 2, 3> Stiffness::gradients =
      Eigen::Matrix<double, 2, 3>
      {
        {-1, 1, 0},
        {-1, 0, 1}
      };

matrix Stiffness::generate_matrix() const {

  auto mesh_ptr = mesh.lock();

  const index_t N_active_vertex = mesh_ptr->n_active_vertices();
  const index_t N_triangles = mesh_ptr->n_triangles();

  //npdebug("Number of active vertices: ", N_active_vertex)
  //npdebug("Number of active triangles: ", N_triangles)
  
  // initialize with a null matrix 
  matrix A = matrix::Zero(N_active_vertex, N_active_vertex);

  for (index_t k = 0; k < N_triangles; ++k) {
    
    // get triangle structure
    const triangle_t& K = mesh.lock()->get_triangle(k);

    // get pre-computed metric structure
    const auto& metric = metrics[k];

    // integrate field
    double field_integral = integrator->integrate(field, metric.transform); 

    //npdebug("Volume: ", metric.abs_detB)
    //npdebug("Field integral: ", field_integral)

    // store triangle contribution
    Eigen::Matrix3d A_k;
    
    // compile triangle matrix
    for(index_t i = 0; i < 3; ++i) {

      // compute metrics
      vertex_t full_grad_i = metric.invBT * gradients.col(i);

      for(index_t j = 0; j < 3; ++j) {

        vertex_t full_grad_j = metric.invBT * gradients.col(j);
        
        // compile contribution
        // int_field(i,j) * < B_k^(-T) * grad(phi_i) , B_k^(-T) * grad(phi_j) > * det(B_k)
        A_k(i, j) = field_integral * metric.abs_detB * full_grad_i.dot(full_grad_j);
      }
    }

    //bool boundary = false;
    //std::vector<index_t> boundary_index;

    // re-iterate and assign to global matrix
    // skip if the vertex is not active
    for (index_t i = 0; i < 3; ++i) {

      // determine target in the edge
      const index_t N_i = K[i];

      // if N_i is out of scope: do nothing
      if (N_i >= N_active_vertex) {
        //boundary_index.push_back(N_i);
        //boundary = true;
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
        A(N_i, N_j) += A_k(i,j);
      }
    }
    
    // debug info
    /*if (boundary) {
      std::cout << "Boundary detected: ";
      for (index_t i : boundary_index)
         std::cout << i << " ";
      std::cout << std::endl;

      npdebug(A_k)
    }*/
  } 
  
  return A;
}

/*
 *  Inherited classes
 */

NeumannStiffness::NeumannStiffness(
          std::weak_ptr<const Mesh> mesh, 
          const ScalarField& field,
          std::shared_ptr<const Integrator> integrator) 
  : Stiffness(std::static_pointer_cast<const AbstractMesh>(mesh.lock()), field, integrator)
{
}


DirichletStiffness::DirichletStiffness(
          std::weak_ptr<const IslandMesh> mesh, 
          const ScalarField& field,
          std::shared_ptr<const Integrator> integrator) 
  : Stiffness(std::static_pointer_cast<const AbstractMesh>(mesh.lock()), field, integrator)
{
}
