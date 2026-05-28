#include "linalg/core/matrix.hpp"

#include <iostream>

// ./build/bin/example_app

int main(int argc, char **argv) {

  std::cout << "Hello world! This is an example executable." << std::endl;

  std::size_t n_rows = 3;
  std::size_t n_cols = 4;

  la::Matrix2D<double> matrix1 = la::zeros(n_rows, n_cols);

  la::Matrix2D<double> matrix2 = la::zeros(n_rows);

  la::Matrix2D<double> matrix3 = la::ones(n_rows, n_cols);

  la::Matrix2D<double> matrix4 = la::ones(n_rows);

  la::Matrix2D<double> matrix5 = la::eye(n_rows);

  matrix1.print();
  matrix2.print();
  matrix3.print();
  matrix4.print();
  matrix5.print();

  // Print the transposed matrices
  matrix1.transpose().print();
  matrix3.transpose().print();

  return 0;
  //   // --- Deprecated examples ---

  //   // ----- Matrix Examples -----
  //   Matrix A = Matrix({{1, 2, 3}, {2, 1, 3}, {-6, -4, 2}});
  //   A.print();

  //   Matrix B = A.inv();

  //   std::cout << "--- B = A^(-1) ---" << std::endl;
  //   B.print();

  //   std::cout << "--- A * B ---" << std::endl;
  //   (A * B).print();

  //   std::cout << "--- B * A ---" << std::endl;
  //   (B * A).print();

  //   std::cout << "--- B^(-1)*v = A*v ---" << std::endl;
  //   util::print(B / util::vec{1, 1, 1});

  //   Matrix Ones = Matrix(util::ones(3));

  //   std::cout << "--- 1_3x3 * 1_3x3 ---" << std::endl;
  //   Ones.mult(Ones).print();

  //   std::cout << "--- A + 1_3x3 ---" << std::endl;
  //   (A + Ones).print();

  //   std::cout << "--- A | 1_3x3 ---" << std::endl;
  //   (A | Ones).print();

  //   std::cout << "--- A + 6*I ---" << std::endl;
  //   (A + (Matrix(util::eye(3)) * 6)).print();

  //   std::cout << "--- (A + 6*I)^(-1) ---" << std::endl;
  //   (A + (Matrix(util::eye(3)) * 6)).inv().print();

  //   // ----- Simplex Example -----
  //   std::cout << "----- Simplex Example -----" << std::endl;

  //   util::matrix A0{{1, 3}, {1, 1}, {3, 1}};
  //   util::vec b0{15, 7, 15};
  //   util::vec c0{2, 2};
  //   // Considering the LP min(x) c^t*x s.t. A*x<=0, x>=0

  //   Simplex model{};

  //   model.addVariables("x", 2, pos);
  //   model.addConstraints("x", A0, b0, leq);
  //   model.addObjective("x", c0, min);

  //   auto [sol, status] = model.solve();

  //   util::print(sol);
}
