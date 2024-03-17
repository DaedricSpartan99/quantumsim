#include "stiffness.hpp"

using namespace qsim2d;

cpx_matrix MassComponent::generate_matrix() const {
  
  // initialize with a null matrix 
  cpx_matrix M = cpx_matrix(N_vertex, point_t(N_vertex), 0);

  for (index_t k = 0; k < N_triangles; ++k) {
    
    for(index_t i = 0; i < interp.size(); ++i) {
      for(index_t j = 0; j < interp.size(); ++j) {
      
        // determine vertex indexes
        const index_t N_i = contributions[k].vert_indexes[i];
        const index_t N_j = contributions[k].vert_indexes[j];

        // compute metrics
        double metric = contributions[k].abs_detB;

        // integrate by interpolation
        const std::complex integral = 0;
      
        for (index_t l = 0; l < interp.size(); ++l) {

          // evaluate local basis
          double phi_il = basis[i](interp[l].location);
          double phi_jl = basis[j](interp[l].location);

          // phi_i * phi_j * det(B_k)
          integral += phi_il * phi_jl * contributions[k].field_evals[l];
        }

        // update matrix
        M[N_i][N_j] += metric * integral;
      }
    }
  } 
  
  return M;
}

