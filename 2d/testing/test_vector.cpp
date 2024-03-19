#include "types.hpp"
#include <iostream>

using namespace qsim2d;

void output_vector(const vertex_t&);

int main() {
  
  using namespace std;

  // define vertex_t constants
  vertex_t a = {1, 1};
  vertex_t b = {2, 3};

  double alpha = 2;

  // vertex_t sum test
  cout << "a + b = ";
  output_vector(a + b);

  // vertex_t scalar multiplication test
  cout << "a * alpha = ";
  output_vector(a * alpha);

  // vertex_t scalar product
  cout << "a * b = ";
  cout << (a * b) << endl;

  // matrix_t sum test
}

void output_vector(const vertex_t& v) {

  using namespace std;
  cout << "(" << v[0] << ", " << v[1] << ")" << endl;
}
