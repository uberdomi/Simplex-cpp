#include "Util.h"
#include <memory>

class Matrix {
    private:
    util::matrix _A; // Immutable
    int _n,_m;

    public:
    enum type{
        eye, ones, zero
    };

    Matrix(util::matrix A);
    Matrix(type, int dim);

    util::vec mult(util::vec v); // A*v
    util::vec multT(util::vec v); // A^t*v
    std::unique_ptr<Matrix> mult(util::matrix B); // A*B
    std::unique_ptr<Matrix> T(); // A^t
    util::vec solve(util::vec v); // A^(-1)*v
    std::unique_ptr<Matrix> inv(); // A^(-1)
    double det(); // det(A)
    util::matrix getMatrix(); // A

    private:
    std::pair<double, util::matrix> gauss(); // perform the gauss algorithm, getting det. and inv. together

};