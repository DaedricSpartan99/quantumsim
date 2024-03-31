#include "matrix_component.hpp"
#include "mesh.hpp"

#include "debug.hpp"

using namespace qsim2d;

MatrixComponent::MatrixComponent(const ScalarField& field, const Mesh& mesh, const std::vector<GaussPair>& points) 
  : N_triangles(mesh.n_triangles()), 
    N_vertex(mesh.n_vertices()), 
    interp(points), 
    contributions(N_triangles) 
{
 
  npdebug("Initializing contributions to size ", interp.size())

  // initialize each triangle and pre-allocate sizes
  for(index_t k = 0; k < N_triangles; ++k) {
    
    contributions[k].field_evals.resize(interp.size());
  }
  
  // inplace mesh and field
  npdebug("Injecting mesh and field")
  update_mesh_and_field(mesh, field); 
}

void MatrixComponent::update_field(const ScalarField& field) {

  npdebug("Is function bound: ", static_cast<bool>(field))

  // update weighted evaluations
  for(index_t k = 0; k < N_triangles; ++k) {

    for (index_t l = 0; l < interp.size(); ++l) {

      npdebug("Update field: ", "compute field eval ", k, " ", l)
      
      // compute global location of interpolation point
      vertex_t location = contributions[k].transform(interp[l].location);

      //npdebug("Location: ", "x = ", location[0], ", y = ", location[1])
      //npdebug("Field: ", field(location))

      // evaluate field in global location and apply interpolation weight
      contributions[k].field_evals[l] = interp[l].weight * field(location);
    }
  }
}

void MatrixComponent::update_mesh_and_field(const Mesh& mesh, const ScalarField& field) {

  N_triangles = mesh.n_triangles();

  // update metrics
  for(index_t k = 0; k < N_triangles; ++k) {

    npdebug("Computing metrics of triangle: ", k)
  
    // get triangle
    triangle_t K = mesh.get_triangle(k);

    const vertex_t a = mesh.get_vertex(K[0]);
    const vertex_t b = mesh.get_vertex(K[1]);
    const vertex_t c = mesh.get_vertex(K[2]);

    // compute matrix (b - a, c - a)
    matrix_t B {
      { b[0] - a[0], c[0] - a[0]}, 
      { b[1] - a[1], c[1] - a[1]}
    };
    
    // define transformation
    contributions[k].transform = [a,B](vertex_t z) -> vertex_t { 

      return a + B * z; 
    };

    // compute determinant
    const double detB = B.determinant();

    // compute inverse transposed of the matrix B
    const matrix_t invBT = B.inverse().transpose();
    
    // store values
    contributions[k].abs_detB = abs(detB);
    contributions[k].invBT = invBT;
    contributions[k].vert_indexes = K;

    npdebug("Volume: ", abs(detB))
  }

  // compute evaluations
  update_field(field); 
}

