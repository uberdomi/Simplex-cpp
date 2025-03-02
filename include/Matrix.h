#ifndef _MATRIX
#define _MATRIX

#include "Utils.h"
#include <memory>

class Matrix {
    private:
    // Immutable
    const util::matrix _A;
    const size_t _n,_m;

    // Can be updated
    util::matrix _A_inv{};
    double _det{0.0};

    public:

    // Constructor

    Matrix(const util::matrix& A);

    // Static constructors

    static Matrix eye(size_t n);
    static Matrix ones(size_t n, size_t m=0);
    static Matrix zeros(size_t n, size_t m=0);

    // Matrix - vector multiplication

    util::vec mult(const util::vec& v) const ; // A*v
    util::vec multT(const util::vec& v) const ; // A^t*v

    Matrix mult(const util::matrix& B) const ; // A*B
    Matrix mult(const Matrix& B) const ; // A*B

    // Solving linear systems

    Matrix T() const ; // A^t
    Matrix inv(); // A^(-1)
    double det(); // det(A)
    util::vec solve(const util::vec& v) const ; // A^(-1)*v

    // Helper functions

    util::matrix getMatrix() const ; // A

    void print() const;

    private:

    std::pair<double, util::matrix> gauss(util::matrix&& RHS) const ; // perform the gauss algorithm, getting det. and inv. together
    std::pair<double, util::matrix> gauss(const util::vec& RHS) const ;

};

#endif // _MATRIX