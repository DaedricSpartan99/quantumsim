#pragma once

#include "types.hpp"
#include "abstract_mesh.hpp"
#include <memory>

namespace qsim2d {

  /*
   * Class grouping all
   */
 
  class FEComponent {

    public:

      FEComponent(std::weak_ptr<const AbstractMesh> mesh,
          const ScalarField& field) : 
            mesh(mesh), 
            metrics(mesh.lock()->compute_metrics()),
            field(field) 
      {
      }


      // Mesh get/set
      inline void set_mesh(std::weak_ptr<const AbstractMesh> mesh) {

        this->mesh = mesh;
        
        // recompute metrics
        metrics = mesh.lock()->compute_metrics();
      }

      inline const std::shared_ptr<const AbstractMesh> get_mesh() const {
        return mesh.lock();
      }

      // Field get/set/eval
      inline void set_field(const ScalarField& field) {
        this->field = field;
      }

      inline ScalarField get_field() const {
        return field;
      }

      inline double eval_field(vertex_t x) const {
        return field(x);
      }
    
    protected:

      // mesh
      std::weak_ptr<const AbstractMesh> mesh;

      // store local contribution
      std::vector<AbstractMesh::triangle_metric> metrics;

      // field
      ScalarField field; 
  };
}
