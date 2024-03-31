#include "mesh.hpp"

using namespace qsim2d;

Mesh::Mesh(const std::vector<vertex_t>& vertices, const std::vector<triangle_t>& triangles) : 
  vertices(vertices), 
  triangles(triangles) 
{
  // TODO: check that it's actuall a mesh
}

Mesh::Mesh(std::vector<vertex_t>&& vertices, std::vector<triangle_t>&& triangles) : 
  vertices(std::move(vertices)), 
  triangles(std::move(triangles)) 
{
  // TODO: check that it's actuall a mesh
}


const vertex_t& Mesh::get_vertex(index_t index) const {
  return vertices[index];
}

vertex_t& Mesh::get_vertex(index_t index) {
  return vertices[index];
}

const std::vector<vertex_t>& Mesh::all_vertices() const {
  return vertices;
}

/*
 * Access to triangles
 */

const triangle_t& Mesh::get_triangle(index_t index) const {
  return triangles[index];
}

triangle_t& Mesh::get_triangle(index_t index) {
  return triangles[index];
}

const std::vector<triangle_t>& Mesh::all_triangles() const {
  return triangles;
}

/*
 * Sizes
 */

index_t Mesh::n_vertices() const {
  return vertices.size();
}

index_t Mesh::n_active_vertices() const {
  return vertices.size();
}


index_t Mesh::n_triangles() const {
  return triangles.size();
}



