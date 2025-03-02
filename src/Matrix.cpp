#include "Matrix.h"
#include <stdexcept>
#include <iostream>

Matrix::Matrix(const util::matrix& A) :
_A(A),
_n{A.size()}, _m{A.at(0).size()}
{
    if(!util::checkRectangular(_A)){
        throw std::invalid_argument("Matrix is not rectangular!");
    }
}

// Static constructors

Matrix Matrix::eye(size_t n) {
    if(n <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }

    util::matrix I(n, util::vec(n, 0.0));

    for(int i=0; i<n; i++) {
        I.at(i).at(i) = 1;
    }

    return Matrix(I);
}

Matrix Matrix::ones(size_t n, size_t m) {
    if(n <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }
    if(m <= 0) {
        m = n;
    }

    util::matrix Ones(n, util::vec(m, 1.0));
    
    return Matrix(Ones);
}

Matrix Matrix::zeros(size_t n, size_t m) {
    if(n <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }
    if(m <= 0) {
        m = n;
    }

    util::matrix Zeros(n, util::vec(m, 0.0));
    
    return Matrix(Zeros);
}

// Matrix - vector multiplication

util::vec Matrix::mult(const util::vec& v) const {
    if(v.size() != _m) {
        throw std::invalid_argument("Matrix and vector have incompatible dimensions!");
    }

    util::vec b(_n, 0.0);

    for(int i=0; i<_n; i++) {
        for(int j=0; j<_m; j++) {
            b.at(i) += _A.at(i).at(j) * v.at(j);
        }
    }

    return b;
}

util::vec Matrix::multT(const util::vec& v) const {
    // Worse cache locality
    if(v.size() != _n) {
        throw std::invalid_argument("Matrix and vector have incompatible dimensions!");
    }

    util::vec b(_m, 0.0);

    for(int j=0; j<_n; j++) {
        for(int i=0; i<_m; i++) {
            b.at(i) += _A.at(j).at(i) * v.at(j);
        }
    }  

    return b;
}

Matrix Matrix::mult(const util::matrix& B) const {
    if(!util::checkRectangular(B)) {
        throw std::invalid_argument("Matrix is not rectangular!");
    }
    if(_m != B.size()) {
        throw std::invalid_argument("Matrices have incompatible dimensions!");
    }

    util::matrix C(_n, util::vec(B.at(0).size(), 0.0));

    for(int i=0; i<_n; i++) {
        for(int j=0; j<_m; j++) {
            for(int k=0; k<B.at(0).size(); k++) {
                C.at(i).at(k) += _A.at(i).at(j) * B.at(j).at(k);
            }
        }
    }
    return Matrix(C);
}

Matrix Matrix::mult(const Matrix& B) const {
    return mult(B.getMatrix());
}

// Solving linear systems

Matrix Matrix::T() const {
    util::matrix B(_m, util::vec(_n, 0.0));

    for(int i=0; i<_n; i++) {
        for(int j=0; j<_m; j++) {
            B.at(j).at(i) = _A.at(i).at(j);
        }
    }

    return Matrix(B);
}

Matrix Matrix::inv() {
    if(_A_inv.empty()) {
        auto [det, Ainv] = gauss(eye(_n).getMatrix());
        _A_inv = Ainv;
        _det = det;
    }
    if(util::isEqual(_det, 0.0)) {
        throw std::runtime_error("Matrix not invertible!");
    }

    return Matrix(_A_inv);
}

double Matrix::det() {
    if(_A_inv.empty()) {
        auto [det, Ainv] = gauss(eye(_n).getMatrix());
        _A_inv = Ainv;
        _det = det;
    }

    return _det;
}

util::vec Matrix::solve(const util::vec& v) const {
    auto [det, Ainv_v] = gauss(v);
    return util::matrixAsVec(Ainv_v);
}

util::matrix Matrix::getMatrix() const {
    return _A;
}

void Matrix::print() const {
    util::print(_A);
}

std::pair<double, util::matrix> Matrix::gauss(util::matrix&& RHS) const {

    // TODO remove printing
    std::cout << "--- Inverting/Solving with matrix: ---" << std::endl;
    util::print(_A);
    std::cout << "-------------------------" << std::endl;

    double det = 1.0;
    util::matrix LHS = _A;
    if(_m != _n) {
        return {0.0, LHS};
    }

    int pivot{0};
    int j{0};
    double tmp{0.0};
    for(int i=0; i<_n; i++) {
        // --- Start Pivot
        // Find biggest entry
        tmp = std::abs(LHS.at(i).at(i));
        pivot=i;
        j=i+1;
        while(j < _n) {
            if(tmp < std::abs(LHS.at(j).at(i))) {
                tmp = std::abs(LHS.at(j).at(i));
                pivot = j;
            }
            j++;
        }

        if(util::isEqual(tmp, 0.0)) {
            // All entries 0, matrix not invertible
            return {0.0, LHS};
        }
        // Swap so that the biggest entry is at the diagonal;
        if(pivot != i) {
            std::swap(LHS.at(i), LHS.at(pivot));
            std::swap(RHS.at(i), RHS.at(pivot));
            det *= (-1);
        }
        // --- End Pivot

        // --- Start Scaling and Subtracting
        // Value to scale, not 0
        tmp = 1 / LHS.at(i).at(i);
        util::scaleRow(LHS.at(i), tmp);
        util::scaleRow(RHS.at(i), tmp);
        det *= tmp;

        j=i+1;
        while(j < _n) {
            tmp = LHS.at(j).at(i);
            // Subtract the scaled i-th row from the j-th row
            util::subtractRow(LHS.at(j), LHS.at(i), tmp);
            util::subtractRow(RHS.at(j), RHS.at(i), tmp);
            j++;
        }
        // --- End Scaling and Subtracting
        
        // TODO remove printing
        std::cout << "------------ LHS -------------" << std::endl;
        util::print(LHS);
        std::cout << "------------ RHS -------------" << std::endl;
        util::print(RHS);
        std::cout << "-------------------------" << std::endl;
    }

    // --- Start Cleanup
    for(int i=_n-1; i>=0; i--) {
        j = i-1;
        while(j >= 0){
            tmp = LHS.at(j).at(i);
            // Subtract the scaled i-th row from the j-th row
            util::subtractRow(LHS.at(j), LHS.at(i), tmp);
            util::subtractRow(RHS.at(j), RHS.at(i), tmp);
            j--;
        }
        // TODO remove printing
        std::cout << "------------ LHS -------------" << std::endl;
        util::print(LHS);
        std::cout << "------------ RHS -------------" << std::endl;
        util::print(RHS);
        std::cout << "-------------------------" << std::endl;
    }
    // --- End Cleanup

    return {det, RHS};
}

std::pair<double, util::matrix> Matrix::gauss(const util::vec& RHS) const {
    return gauss(util::vecAsMatrix(RHS));
}