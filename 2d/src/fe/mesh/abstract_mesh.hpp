#pragma once

#include "types.hpp"
#include <vector>

namespace qsim2d {
  
  /*
   *  Layer of abstraction for a mesh with multiple attached functionalities
   *
   *  This object is conceived to own the active vertices at the beginning:
   *
   *    [active_vertices, inactive_vertices]
   *
   *  By active vertex is meant that it doesn't contribute to the solution of the problem.
   *  For example in the Dirichlet boundary conditions, active vertices are the internal ones. 
   */  
  class AbstractMesh {

    public:
    
      /*
       *  Access to vertices
       */

      virtual const vertex_t& get_vertex(index_t index) const = 0;

      virtual vertex_t& get_vertex(index_t index) = 0;

      /*
       * Access to triangles
       */

      virtual const triangle_t& get_triangle(index_t index) const = 0;

      virtual triangle_t& get_triangle(index_t index) = 0;
      
      /*
       * Size of the mesh
       */      
      
      virtual index_t n_vertices() const = 0;
      virtual index_t n_active_vertices() const = 0;
      virtual index_t n_triangles() const = 0;

      /*
       * Metrics computation for extended classes
       */

      struct triangle_metric {
            
        // transform from local space to global
        // T_k(z) = a + B_k * z
        transform_t transform;
     
        // determinant of B_k
        double abs_detB;
        
        // inverse of B_k for gradient transformation
        matrix_t invBT;
      };

      std::vector<triangle_metric> compute_metrics() const;
  };
}
