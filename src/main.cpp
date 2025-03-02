#include "Utils.h"
#include "Matrix.h"
#include <iostream>

int main(int argc, char** argv) {
    // util::matrix A{{1,2,3}, {2,1,3}, {7,4,2}};
    // util::print(A);

    Matrix A = Matrix({{1,2,3}, {2,1,3}, {-6,-4,2}});
    A.print();

    Matrix B = A.inv();

    std::cout << "--- B = A^(-1) ---" << std::endl;
    B.print();

    std::cout << "--- A * B ---" << std::endl;
    (A * B).print();

    std::cout << "--- B * A ---" << std::endl;
    (B * A).print();

    std::cout << "--- B^(-1)*v = A*v ---" << std::endl;
    util::print(B / util::vec{1,1,1});

    size_t n{4};
    Matrix Ones = Matrix::ones(n);

    Ones.mult(Ones).print();

    // Ones.inv().print();

    return 0;
}