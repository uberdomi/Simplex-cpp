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

    Matrix(const util::matrix& A);
    Matrix(InitType type, size_t dim);

    util::vec mult(const util::vec& v) const ; // A*v
    util::vec multT(const util::vec& v) const ; // A^t*v

    std::unique_ptr<Matrix> mult(const util::matrix& B) const ; // A*B
    std::unique_ptr<Matrix> mult(const Matrix& B) const ; // A*B

    std::unique_ptr<Matrix> T() const ; // A^t
    util::vec solve(const util::vec& v) const ; // A^(-1)*v
    std::unique_ptr<Matrix> inv() const ; // A^(-1)
    double det() const ; // det(A)
    util::matrix getMatrix() const ; // A

    void print() const;

    private:
    static util::matrix initMatrix(InitType type, size_t dim);

    std::pair<double, util::matrix> gauss() const ; // perform the gauss algorithm, getting det. and inv. together

};

#endif // _MATRIX