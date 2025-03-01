#include "Util.h"
#include "Matrix.h"

int main(int argc, char** argv) {
    util::matrix A{{1,2,3}, {2,1,3}, {7,4,2}};
    util::print(A);

    size_t n{4};
    Matrix Ones(Matrix::InitType::ones, n);

    Ones.mult(Ones)->print();
    return 0;
}