#include "load.hpp"

#include "debug.hpp"

using namespace qsim2d;

std::array<std::function<double(const vertex_t&)>, 3> LoadComponent::basis = {
        [](const vertex_t& z) -> double { return 1 - z[0] - z[1]; }, 
        [](const vertex_t& z) -> double { return z[0]; },
        [](const vertex_t& z) -> double { return z[1]; }
      };

cpx_matrix LoadComponent::generate_matrix() const {
  
  // initialize with a null vector
  cpx_matrix M = cpx_matrix::Zero(N_vertex, 1);

  for (index_t k = 0; k < N_triangles; ++k) {
    
    for(index_t i = 0; i < 3; ++i) {
      
        // determine vertex indexes
        const index_t N_i = contributions[k].vert_indexes[i];

        // compute metrics
        double metric = contributions[k].abs_detB;

        // integrate by interpolation
        qsim2d::complex integral = 0.;
      
        for (index_t l = 0; l < interp.size(); ++l) {

          // evaluate local basis
          double phi_il = basis[i](interp[l].location);

          // phi_i * det(B_k)
          integral += phi_il * contributions[k].field_evals[l];
        }

        npdebug("Integral field phi_", i, ": ", integral)

        // update matrix
        M(N_i, 0) += metric * integral;
    }
  } 
  
  return M;
}

