#ifndef _MATRIX
#define _MATRIX

#include "Utils.h"
#include <memory>

class Matrix {
    private:
    // Immutable - stupid, let's make it changeable
    util::matrix _A;
    size_t _n,_m;

    // Can be updated
    util::matrix _A_inv{};
    double _det{0.0};

    // Construct matrix with its inverse
    Matrix(util::matrix&& A, const util::matrix& A_inv, double&& det);

    public:

    // Constructor

    Matrix(const util::matrix& A);

    // Matrix addition

    Matrix add(const util::matrix& B) const ; // A+B
    Matrix add(const Matrix& B) const ; // A+B

    Matrix operator+(const util::matrix& B) const ; // A+B
    Matrix operator+(const Matrix& B) const ; // A+B

    // Matrix concatenation

    Matrix conc(const util::matrix& B) const ; // [A | B]
    Matrix conc(const Matrix& B) const ; // [A | B]

    Matrix operator|(const util::matrix& B) const ; // [A | B]
    Matrix operator|(const Matrix& B) const ; // [A | B]

    // Matrix scaling

    Matrix scale(const double& sf) const ; // sf * A
    Matrix operator*(const double& sf) const ; // sf * A

    // Matrix - Vector multiplication

    util::vec mult(const util::vec& v) const ; // A*v
    util::vec multT(const util::vec& v) const ; // A^t*v

    // Matrix - Matrix multiplication

    Matrix mult(const util::matrix& B) const ; // A*B
    Matrix mult(const Matrix& B) const ; // A*B

    util::vec operator*(const util::vec& v) const ; // A*v
    Matrix operator*(const util::matrix& B) const ; // A*B
    Matrix operator*(const Matrix& B) const ; // A*B

    // Solving linear systems

    Matrix T() const ; // A^t
    Matrix inv(); // A^(-1)
    double det(); // det(A)
    util::vec solve(const util::vec& v) const ; // A^(-1)*v
    util::vec operator/(const util::vec& v) const ; // A^(-1)*v

    Matrix ShermanMorrison(const util::vec& u, const util::vec& v);

    // Helper functions

    util::matrix getMatrix() const ; // A
    util::vec getRow(const int& i) const ; // A_i

    static void swapRows(Matrix& A, Matrix& B, const int& i, const int& j);

    void swapToBack(const int& i);

    void print() const;

    private:

    std::pair<double, util::matrix> gauss(util::matrix&& RHS) const ; // perform the gauss algorithm, getting det. and inv. together
    std::pair<double, util::matrix> gauss(const util::vec& RHS) const ;

};

#endif // _MATRIX