#ifndef _MATRIX
#define _MATRIX

#include "Util.h"
#include <memory>

class Matrix {
    private:
    // Immutable
    const util::matrix _A;
    const size_t _n,_m;

    public:
    enum InitType{
        eye, ones, zeros
    };

    Matrix(util::matrix A);
    Matrix(InitType type, size_t dim);

    util::vec mult(util::vec v); // A*v
    util::vec multT(util::vec v); // A^t*v
    std::unique_ptr<Matrix> mult(util::matrix B); // A*B
    std::unique_ptr<Matrix> T(); // A^t
    util::vec solve(util::vec v); // A^(-1)*v
    std::unique_ptr<Matrix> inv(); // A^(-1)
    double det(); // det(A)
    util::matrix getMatrix(); // A

    private:
    static util::matrix initMatrix(InitType type, size_t dim);

    std::pair<double, util::matrix> gauss(); // perform the gauss algorithm, getting det. and inv. together

};

#endif // _MATRIX