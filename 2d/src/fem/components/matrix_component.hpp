#pragma once

#include "types.hpp"
#include "interpolation.hpp"

namespace qsim2d {

  class Mesh;
  class ScalarField;

  /*
   * Class grouping all information needed to solve stiffness parameter
   * of a Poisson-like equation
   */
  
  class MatrixComponent {

    public:
      MatrixComponent(const ScalarField&, const Mesh&, const std::vector<GaussPair>& gauss = GAUSS_TRI_POINTS);

      virtual cpx_matrix generate_matrix() const = 0;
      
      // update mesh structure: slow
      virtual void update_mesh_and_field(const Mesh&, const ScalarField&);

      // update input field: fast
      virtual void update_field(const ScalarField&);


    protected:
      
      struct triangle_contribution {
            
        // transform from local space to global
        // T_k(z) = a + B_k * z
        transform_t transform;
     
        // determinant of B_k
        double abs_detB;
        
        // inverse of B_k
        matrix_t invBT;
     
        // store vertex indexes
        triangle_t vert_indexes;
     
        // field evaluations for each gauss point
        // weight_l * field(T_k(z_l))
        cpx_vector field_evals;
      };

      
      // triangles size
      std::size_t N_triangles;

      // number of vertices
      index_t N_vertex;

      // store local contribution
      std::vector<triangle_contribution> contributions;
      
      // Local interpolation points
      std::vector<GaussPair> interp;
  };
}
