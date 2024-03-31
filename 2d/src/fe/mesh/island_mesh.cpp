#include "island_mesh.hpp"
#include "mesh.hpp"

#include <unordered_map>
#include <memory>

#include "debug.hpp"

#include <boost/functional/hash.hpp>

using namespace qsim2d;


template<class T>
struct equal_undirected {

  constexpr bool operator()(const T& lhs, const T& rhs) const {
    return 
      (lhs.first == rhs.first and lhs.second == rhs.second) or
      (lhs.first == rhs.second and lhs.second == rhs.first);
  }
};


// Only for pairs of std::hash-able types for simplicity.
// You can of course template this struct to allow other hash functions
template<class Pair>
struct hash_pair {
   std::size_t operator()(const Pair &p) const noexcept {
      std::size_t seed = 0;
      boost::hash_combine(seed, p.first);
      boost::hash_combine(seed, p.second);
      return seed;
    }
};

std::vector<bool> separe_internal_and_boundary(const std::vector<triangle_t>& triangles, index_t Nv) {

  /*
   * Identify internal and boundary vertices
   *
   * Algorithm:
   *  - Construct edges from triangles
   *  - For each edge: count the number of polygons that own the edge
   *  - Take only edges that have count == 1
   *  - Extract vertices set
   *  - Triangles owning those edges are marked as boundary triangles
   */

  typedef std::pair<index_t, index_t> edge;

  // collect edges
  std::unordered_map<edge, int, hash_pair<edge>, equal_undirected<edge>> edge_counter;
  
  for (size_t k = 0; k < triangles.size(); ++k) {

    const triangle_t& T = triangles[k];
    
    // construct and insert edges
    edge e_01(T[0], T[1]);
    edge e_12(T[1], T[2]);
    edge e_20(T[2], T[0]);
  
    edge_counter[e_01]++;
    edge_counter[e_12]++;
    edge_counter[e_20]++;
  }

  // detect internal edges
  std::vector<edge> boundary_edges;

  for (auto const& [key, count] : edge_counter) {
    if (count <= 1) {
      boundary_edges.push_back(key);
    }
  }

  // assign internal vertices
  std::vector<bool> is_internal(Nv, true);
  
  for (edge const& e : boundary_edges) {
    is_internal[e.first] = false;
    is_internal[e.second] = false;
  } 

  return is_internal;
}

IslandMesh::IslandMesh(const std::vector<vertex_t>& vertices, const std::vector<triangle_t>& triangles) :
  internal_mesh(nullptr),
  boundary_vertices(),  // init empty
  boundary_triangles()  // init empty
{
  auto Nv = vertices.size();
  auto Nk = triangles.size();
  
  // get internal vertices 
  std::vector<bool> is_internal(separe_internal_and_boundary(triangles, Nv));

  // initialize internal components 
  std::vector<vertex_t> internal_vertices;
  std::vector<triangle_t> internal_triangles;
  
  // store pivoting information
  std::vector<int> vertex_pivot(Nv, -1);

  // re-order vertices
  for (index_t i = 0; i < Nv; ++i) {
    
    if (not is_internal[i]) {
      
      // boundary

      // store pivoting
      vertex_pivot[i] = boundary_vertices.size();

      // add boundary, already marked
      boundary_vertices.push_back(vertices[i]);

    } else {

      // internal

      // store pivoting
      vertex_pivot[i] = internal_vertices.size();

      // add to internal
      internal_vertices.push_back(vertices[i]);
    }

    //npdebug("Vertex ", i, " remapped to ", vertex_pivot[i], ". Internal: ", is_internal[i])
  }
  
  /*
   *  Construct sub-mesh for boundary conditions.
   *
   *  Remap vertex indexing in order to match the structure:
   *
   *  [internal_nodes, boundary_nodes]
   *
   *  Triangles must be re-mapped in order to match the new vertex intexing.
   *  Some boundary triangles will point to internal vertices.
   */
  const size_t Nvi = internal_vertices.size();

  // use marked vertices to identify triangles
  for (size_t k = 0; k < Nk; ++k) {
    
    // check if all points are internal
    bool internal = true;

    for (index_t j = 0; j < 3; ++j) {
      internal &= is_internal[triangles[k][j]];
    }
    
    // remap triangles with new indexing 
    triangle_t remapped_triangle;

    for (index_t j = 0; j < 3; ++j) { 

      // if on boundary (vertexes), shift after the end
      index_t shift = (is_internal[triangles[k][j]]) ? 0 : Nvi;

      // remap
      remapped_triangle[j] = vertex_pivot[triangles[k][j]] + shift;
    }
    
    // debug info
    /*npdebug("Triangle ", 
        triangles[k][0], " ", triangles[k][1], " ", triangles[k][2], " remapped to ", 
        remapped_triangle[0], " ", remapped_triangle[1], " ", remapped_triangle[2], ", Internal: ", internal
        )*/

    // if it's the case, mark triangle as internal
    if (internal) {
      internal_triangles.push_back(remapped_triangle);
    } else {
      boundary_triangles.push_back(remapped_triangle);
    }
  }

  // finally construct internal mesh
  internal_mesh = std::make_unique<Mesh>(std::move(internal_vertices), std::move(internal_triangles));
}
    
 /*
  *  Access to vertices
  */

const vertex_t& IslandMesh::get_vertex(index_t index) const {

  return (index < internal_mesh->n_vertices()) ? internal_mesh->get_vertex(index) : boundary_vertices[index - internal_mesh->n_vertices()];
}

vertex_t& IslandMesh::get_vertex(index_t index) {

  return (index < internal_mesh->n_vertices()) ? internal_mesh->get_vertex(index) : boundary_vertices[index - internal_mesh->n_vertices()];
}

 /*
  * Access to triangles
  */

const triangle_t& IslandMesh::get_triangle(index_t index) const {

  return (index < internal_mesh->n_triangles()) ? internal_mesh->get_triangle(index) : boundary_triangles[index - internal_mesh->n_triangles()];
}

triangle_t& IslandMesh::get_triangle(index_t index) {

  return (index < internal_mesh->n_triangles()) ? internal_mesh->get_triangle(index) : boundary_triangles[index - internal_mesh->n_triangles()];
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

index_t IslandMesh::n_active_vertices() const {
  return internal_mesh->n_vertices();
}

index_t IslandMesh::n_triangles() const {
  return boundary_triangles.size() + internal_mesh->n_triangles();
}



