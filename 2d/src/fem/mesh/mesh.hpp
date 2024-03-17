#pragma once

#include "abstract_mesh.hpp"

namespace qsim2d {
  
  class Mesh : public AbstractMesh {

    public:
    
      Mesh(const std::vector<vertex_t>&, const std::vector<triangle_t>&);
      Mesh(std::vector<vertex_t>&&, std::vector<triangle_t>&&);
    
      /*
       *  Access to vertices
       */

      virtual const vertex_t& get_vertex(index_t index) const override;

      virtual vertex_t& get_vertex(index_t index) override;

      const std::vector<vertex_t>& all_vertices() const;

      /*
       * Access to triangles
       */

      virtual const triangle_t& get_triangle(index_t index) const override;

      virtual triangle_t& get_triangle(index_t index) override;
      
      const std::vector<vertex_t>& all_triangles() const;

      /*
       * Sizes
       */

      index_t n_vertices() const override;
      index_t n_triangles() const override;

    private:

        std::vector<vertex_t> vertices;
        std::vector<triangle_t> triangles;
  };
}
