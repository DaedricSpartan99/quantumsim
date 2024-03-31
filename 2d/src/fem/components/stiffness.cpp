#include "stiffness.hpp"

#include "debug.hpp"

using namespace qsim2d;

const Eigen::Matrix<double, 2, 3> StiffnessComponent::gradients
{
  {-1, 1, 0},
  {-1, 0, 1}
};

cpx_matrix StiffnessComponent::generate_matrix() const {
  
  // initialize with a null matrix 
  cpx_matrix A = cpx_matrix::Zero(N_vertex, N_vertex);

  for (index_t k = 0; k < N_triangles; ++k) {
    
    for(index_t i = 0; i < 3; ++i) {
      for(index_t j = 0; j < 3; ++j) {
      
        // determine vertex indexes
        const index_t N_i = contributions[k].vert_indexes[i];
        const index_t N_j = contributions[k].vert_indexes[j];

        // compute metrics
        // < B_k^(-T) * grad(phi_i) , B_k^(-T) * grad(phi_j) > * det(B_k)
        vertex_t full_grad_i = contributions[k].invBT * gradients.col(i);
        vertex_t full_grad_j = contributions[k].invBT * gradients.col(j);
        double metric = contributions[k].abs_detB * full_grad_i.dot(full_grad_j);

        // integrate by interpolation
        qsim2d::complex integral = 0;
      
        for (index_t l = 0; l < interp.size(); ++l) {
          integral += contributions[k].field_evals[l];
        }

        npdebug("Integral stiffness: ", integral.real())
        npdebug("Metric stiffness: ", metric)

        // update matrix
        A(N_i, N_j) += metric * integral;
      }
    }
  } 
  
  return A;
}


