#pragma once

#include "types.hpp"
#include <vector>

namespace qsim2d {
  
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
      virtual index_t n_triangles() const = 0;
  };
}
