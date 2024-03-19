#include "mass.hpp"

using namespace qsim2d;

std::array<std::function<double(const vertex_t&)>, 3> MassComponent::basis = {
        [](const vertex_t& z) -> double { return 1 - z[0] - z[1]; }, 
        [](const vertex_t& z) -> double { return z[0]; },
        [](const vertex_t& z) -> double { return z[1]; }
      };

cpx_matrix MassComponent::generate_matrix() const {
  
  // initialize with a null matrix 
  cpx_matrix M = cpx_matrix(N_vertex, cpx_vector(N_vertex, 0));

  for (index_t k = 0; k < N_triangles; ++k) {
    
    for(index_t i = 0; i < interp.size(); ++i) {
      for(index_t j = 0; j < interp.size(); ++j) {
      
        // determine vertex indexes
        const index_t N_i = contributions[k].vert_indexes[i];
        const index_t N_j = contributions[k].vert_indexes[j];

        // compute metrics
        double metric = contributions[k].abs_detB;

        // integrate by interpolation
        qsim2d::complex integral = 0.;
      
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

