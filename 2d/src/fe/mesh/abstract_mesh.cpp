#include "abstract_mesh.hpp"

#include "debug.hpp"

using namespace qsim2d;

std::vector<AbstractMesh::triangle_metric> AbstractMesh::compute_metrics() const {
  
  const index_t N_triangles = n_triangles();

  std::vector<AbstractMesh::triangle_metric> contributions;
  contributions.reserve(N_triangles);

  for(index_t k = 0; k < N_triangles; ++k) {

    // get triangle
    const triangle_t& K = get_triangle(k);

    const vertex_t& a = get_vertex(K[0]);
    const vertex_t& b = get_vertex(K[1]);
    const vertex_t& c = get_vertex(K[2]);

    
    // compute matrix (b - a, c - a)
    matrix_t B {
      { b[0] - a[0], c[0] - a[0]}, 
      { b[1] - a[1], c[1] - a[1]}
    };

    triangle_metric metric;
    
    // define transformation
    metric.transform = [a,B](vertex_t z) -> vertex_t { 

      return a + B * z; 
    };

    // compute determinant
    const double detB = B.determinant();
    
    // debug info
    /*if (abs(detB) < 1e-15 ) {
      npdebug("Volume zero detected: ")
      npdebug("Triangle: ", K[0], " ", K[1], " ", K[2])
      npdebug("a = ", a)
      npdebug("b = ", b)
      npdebug("c = ", c)
      npdebug()
    }*/

    // compute inverse transposed of the matrix B
    const matrix_t invBT = B.inverse().transpose();
    
    // store values
    metric.abs_detB = abs(detB);
    metric.invBT = invBT;
   
    contributions.push_back(metric); 
  }

  return contributions;
}
