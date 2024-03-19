#include "stiffness.hpp"

using namespace qsim2d;

cpx_matrix StiffnessComponent::generate_matrix() const {
  
  // initialize with a null matrix 
  cpx_matrix A = cpx_matrix(N_vertex, cpx_vector(N_vertex, 0));

  for (index_t k = 0; k < N_triangles; ++k) {
    
    for(index_t i = 0; i < interp.size(); ++i) {
      for(index_t j = 0; j < interp.size(); ++j) {
      
        // determine vertex indexes
        const index_t N_i = contributions[k].vert_indexes[i];
        const index_t N_j = contributions[k].vert_indexes[j];

        // compute metrics
        // < B_k^(-T) * grad(phi_i) , B_k^(-T) * grad(phi_j) > * det(B_k)
        vertex_t full_grad_i = contributions[k].invBT * gradients[i];
        vertex_t full_grad_j = contributions[k].invBT * gradients[j];
        double metric = contributions[k].abs_detB * (full_grad_i * full_grad_j);

        // integrate by interpolation
        qsim2d::complex integral = 0;
      
        for (index_t l = 0; l < interp.size(); ++l) {
          integral += contributions[k].field_evals[l];
        }

        // update matrix
        A[N_i][N_j] += metric * integral;
      }
    }
  } 
  
  return A;
}


