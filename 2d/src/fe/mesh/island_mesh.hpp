#pragma once

#include <memory>
#include "abstract_mesh.hpp"

namespace qsim2d {

  class Mesh;

  /*
   * A mesh that automatically separes internal nodes from boundary nodes.
   * 
   * The name "island" refers to a bounded connected domain, where the boundary is the "beach" and the 
   * internal part is the "land".
   */
  
  class IslandMesh : public AbstractMesh {

    public:
    
      IslandMesh(const std::vector<vertex_t>&, const std::vector<triangle_t>&);
    
      /*
       *  Access to vertices
       *
       *  Careful: can be unefficient due to cache misses if jumping from internal to external vertices
       */

      virtual const vertex_t& get_vertex(index_t index) const override;

      virtual vertex_t& get_vertex(index_t index) override;

      /*
       * Access to triangles
       *
       *  Careful: can be unefficient due to cache misses if jumping from internal to external triangles
       */

      virtual const triangle_t& get_triangle(index_t index) const override;

      virtual triangle_t& get_triangle(index_t index) override;

      /*
       * Const access to components
       */
      
      const Mesh& get_internal_mesh() const;  

      const std::vector<vertex_t>& get_boundary_vertices() const;
      const std::vector<triangle_t>& get_boundary_triangles() const;

      /*
       * Sizes
       */

      index_t n_vertices() const override;
      index_t n_active_vertices() const override;
      index_t n_triangles() const override;

    private:
        
      // store internal values
      std::unique_ptr<Mesh> internal_mesh;
      
      // store boundary information
      std::vector<vertex_t> boundary_vertices;
      std::vector<triangle_t> boundary_triangles;
  };
}
