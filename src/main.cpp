#include "Utils.h"
#include "Matrix.h"

int main(int argc, char** argv) {
    // util::matrix A{{1,2,3}, {2,1,3}, {7,4,2}};
    // util::print(A);

    Matrix A = Matrix({{1,2,3}, {2,1,3}, {-6,-4,2}});
    A.print();

    Matrix B = A.inv();

    B.print();

    A.mult(B).print();

    B.mult(A).print();

    size_t n{4};
    Matrix Ones = Matrix::ones(n);

    Ones.mult(Ones).print();

    // Ones.inv().print();

    return 0;
}