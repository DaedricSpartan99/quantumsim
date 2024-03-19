#include "matrix_component.hpp"
#include "mesh.hpp"
#include "scalar_field.hpp"

using namespace qsim2d;

MatrixComponent::MatrixComponent(const ScalarField& field, const Mesh& mesh, const std::vector<GaussPair>& points) 
  : N_triangles(mesh.n_triangles()), 
    N_vertex(mesh.n_vertices()), 
    interp(points), 
    contributions(N_triangles) 
{
 
  // initialize each triangle and pre-allocate sizes
  for(index_t k = 0; k < N_triangles; ++k) {

    contributions[k].field_evals.reserve(interp.size());
  }
  
  // inplace mesh and field
  update_mesh_and_field(mesh, field); 
}

void MatrixComponent::update_field(const ScalarField& field) {

  // update weighted evaluations
  for(index_t k = 0; k < N_triangles; ++k) {

    for (index_t l = 0; l < interp.size(); ++l) {
      
      // compute global location of interpolation point
      auto location = contributions[k].transform(interp[l].location);

      // evaluate field in global location and apply interpolation weight
      contributions[k].field_evals[l] = interp[l].weight * field(location);
    }
  }
}

void MatrixComponent::update_mesh_and_field(const Mesh& mesh, const ScalarField& field) {

  N_triangles = mesh.n_triangles();

  // update metrics
  for(index_t k = 0; k < N_triangles; ++k) {
  
    // get triangle
    triangle_t K = mesh.get_triangle(k);

    const vertex_t a = mesh.get_vertex(K[0]);
    const vertex_t b = mesh.get_vertex(K[1]);
    const vertex_t c = mesh.get_vertex(K[2]);

    // compute matrix (b - a, c - a)
    matrix_t B = {{
      { b[0] - a[0], c[0] - a[0]}, 
      { b[1] - a[1], c[1] - a[1]}
    }};
    
    // define transformation
    contributions[k].transform = [a,B](vertex_t z) -> vertex_t { 

      return a + B * z; 
    };

    // compute determinant
    const double detB = B[0][0] * B[1][1] - B[1][0] * B[0][1];

    // compute inverse transposed of the matrix B
    const matrix_t invBT = {{
      {B[1][1] / detB, -B[1][0] / detB}, 
      {-B[0][1] / detB , B[0][0] / detB}
    }};
    
    // store values
    contributions[k].abs_detB = abs(detB);
    contributions[k].invBT = invBT;
    contributions[k].vert_indexes = K;
  }

  // compute evaluations
  update_field(field); 
}

