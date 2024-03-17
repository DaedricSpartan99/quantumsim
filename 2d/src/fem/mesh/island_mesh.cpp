#include "island_mesh.hpp"
#include "mesh.hpp"

IslandMesh::IslandMesh(const std::vector<vertex_t>& vertices, const std::vector<triangle_t>& triangles) :
  internal_mesh(nullptr),
  boundary_vertices(),  // init empty
  boundary_triangles()  // init empty
{
  auto Nv = vertices.size();
  auto Nk = triangles.size();
  
  // loop over triangles and store occurrence count for vertices
  std::vector<int> counts(Nv, 0);

  for (size_t k = 0; k < Nk; ++k) {

    for (index_t j = 0; j < 3; ++j)
      // increase counts
      ++counts[triangles[k][j]];
  }
  
  // initialize internal components 
  std::vector<bool> is_internal(Nv, false);
  std::vector<vertex_t> internal_vertices;
  std::vector<triangle_t> internal_triangles;

  // identify vertices that have 2 or less counts
  // and set them to boundary, otherwise they're internal
  for (index_t i = 0; i < Nv; ++i) {

    if (counts[i] <= 2) {

      // add boundary, already marked
      boundary_vertices.push(vertices[i]);

    } else {

      // add to internal
      internal_vertices.push(vertices[i]);

      // mark as internal
      is_internal[i] = true;
    }
  }

  // use marked vertices to identify triangles
  for (size_t k = 0; k < Nk; ++k) {
    
    // check if whether any internal points is present
    bool internal = false;

    for (index_t j = 0; j < 3; ++j) {
      internal |= is_internal[triangles[k][j]];
    }
    
    // if it's the case, mark triangle as internal
    if (internal) 
      internal_triangles.push_back(triangles[k]);
    else
      boundary_triangles.push_back(triangles[k]);
  }

  // finally construct internal mesh
  internal_mesh = std::make_unique<Mesh>(std::move(internal_vertices), std::move(boundary_vertices));
}
    
 /*
  *  Access to vertices
  */

const vertex_t& IslandMesh::get_vertex(index_t index) {

  return (index < internal_mesh->n_vertices()) mesh.get_vertex(index) : boundary_vertices[index - internal_mesh->n_vertices()];
}

vertex_t& IslandMesh::get_vertex(index_t index) {

  return (index < internal_mesh->n_vertices()) mesh.get_vertex(index) : boundary_vertices[index - internal_mesh->n_vertices()];
}

 /*
  * Access to triangles
  */

const triangle_t& IslandMesh::get_triangle(index_t index) const {

  return (index < internal_mesh->n_triangles()) mesh.get_triangle(index) : boundary_triangles[index - internal_mesh->n_triangles()];
}

triangle_t& IslandMesh::get_triangle(index_t index) {

  return (index < internal_mesh->n_triangles()) mesh.get_triangle(index) : boundary_triangles[index - internal_mesh->n_triangles()];
}

 /*
  * Const access to components
  */
 
const Mesh& IslandMesh::get_internal_mesh() const {
  return *internal_mesh;
}

const std::vector<vertex_t>& IslandMesh::get_boundary_vertices() const {
  return boundary_vertices;
}

const std::vector<triangle_t>& IslandMesh::get_boundary_triangles() const {
  return boundary_triangles;
}

/*
 * Sizes
 */

index_t IslandMesh::n_vertices() const {
  return boundary_vertices.size() + internal_mesh->n_vertices();
}

index_t IslandMesh::n_triangles() const {
  return boundary_triangles.size() + internal_mesh->n_triangles();
}



