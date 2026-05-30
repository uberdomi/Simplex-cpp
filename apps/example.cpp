#include "linalg/core/matrix.hpp"
#include "linalg/ops/basic.hpp"

#include <iostream>
#include <string>
#include <vector>

// Compile
// cmake --build build --target example_app --
// Run
// ./build/bin/example_app

int main(int argc, char **argv) {

  std::cout << "Hello world! This is an example executable." << std::endl;

  // --- Matrix instantiations ---

  std::size_t n_rows = 3;
  std::size_t n_cols = 4;

  la::Matrix2D matrix1 = la::zeros(n_rows, n_cols);

  la::Matrix2D matrix2 = la::zeros(n_rows);

  la::Matrix2D matrix3 = la::ones(n_rows, n_cols);

  la::Matrix2D matrix4 = la::ones(n_rows);

  la::Matrix2D matrix5 = la::eye(n_rows);

  matrix1.print();
  matrix2.print();
  matrix3.print();
  matrix4.print();
  matrix5.print();

  // Print the transposed matrices
  std::cout << "Transposed matrices and their presumed shapes: " << std::endl;
  auto t1 = matrix1.transpose();
  t1.print();
  std::cout << "T1 : (" << std::to_string(t1.get_n_rows()) << ","
            << std::to_string(t1.get_n_cols()) << ")" << std::endl;

  auto t3 = matrix1.transpose();
  t3.print();
  std::cout << "T3 : (" << std::to_string(t3.get_n_rows()) << ","
            << std::to_string(t3.get_n_cols()) << ")" << std::endl;

  // Check the originals again
  matrix1.print();
  std::cout << "M1 : (" << std::to_string(matrix1.get_n_rows()) << ","
            << std::to_string(matrix1.get_n_cols()) << ")" << std::endl;

  matrix3.print();
  std::cout << "M3 : (" << std::to_string(matrix3.get_n_rows()) << ","
            << std::to_string(matrix3.get_n_cols()) << ")" << std::endl;

  // --- Basic operations ---
  auto A1 = la::Matrix2D(
      std::vector<std::vector<double>>{{1, 2, 3}, {4, 5, 6}, {6, 6, 7}});

  auto A2 = la::Matrix2D(
      std::vector<std::vector<double>>{{2, 1, 3}, {7, 6, 9}, {0, 0, 1}});

  auto v1 = std::vector<double>{1, 2, 3};
  auto v2 = std::vector<double>{1, 1, 2};

  std::cout << "--- v1 ---" << std::endl;
  la::print(v1);
  std::cout << "--- v2 ---" << std::endl;
  la::print(v2);

  // Vector operations
  std::cout << "--- v1 + v2 ---" << std::endl;
  la::print(la::add(v1, v2));

  std::cout << "--- v1 - v2 ---" << std::endl;
  la::print(la::sub(v1, v2));

  std::cout << "--- v1 * v2 ---" << std::endl;
  la::print(la::mult(v1, v2));

  std::cout << "--- <v1, v2> ---" << std::endl;
  std::cout << std::to_string(la::dot(v1, v2)) << std::endl;

  std::cout << "--- outer(v1, v2) ---" << std::endl;
  la::dyadic(v1, v2).print();

  // Matrix operations
  std::cout << "--- A1 ---" << std::endl;
  A1.print();
  std::cout << "--- A2 ---" << std::endl;
  A2.print();

  std::cout << "--- sum(A1) = " + std::to_string(la::sum(A1)) + " ---"
            << std::endl;

  std::cout << "--- sum(A2) = " + std::to_string(la::sum(A2)) + " ---"
            << std::endl;

  std::cout << "--- row_sum(A1) ---" << std::endl;
  la::print(la::row_sum(A1));

  std::cout << "--- row_sum(A2) ---" << std::endl;
  la::print(la::row_sum(A2));

  std::cout << "--- col_sum(A1) ---" << std::endl;
  la::print(la::col_sum(A1));

  std::cout << "--- col_sum(A2) ---" << std::endl;
  la::print(la::col_sum(A2));

  std::cout << "--- A1 * v1 ---" << std::endl;
  la::print(A1 * v1);

  std::cout << "--- A2 * v2 ---" << std::endl;
  la::print(A2 * v2);

  std::cout << "--- A1 + A2 ---" << std::endl;
  (A1 + A2).print();

  std::cout << "--- A1 - A2 ---" << std::endl;
  (A1 - A2).print();

  std::cout << "--- A1 * A2 ---" << std::endl;
  (A1 * A2).print();

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
