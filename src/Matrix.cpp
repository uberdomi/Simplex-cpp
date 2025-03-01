#include "Matrix.h"
#include <stdexcept>

util::matrix Matrix::initMatrix(InitType type, size_t dim) {
    switch(type) {
        case eye: {
            util::matrix I(dim, util::vec(dim, 0.0));
            for(int i=0; i<dim; i++) {
                I.at(i).at(i) = 1;
            }
            return I;
        };
        case ones: {
            return util::matrix(dim, util::vec(dim, 1.0));
        };
        case zeros: {
            return util::matrix(dim, util::vec(dim, 0.0));
        };
        default : {
            throw std::invalid_argument("Cannot instantiate matrix!");
        }
    }
}

Matrix::Matrix(const util::matrix& A) :
_A(A),
_n{A.size()}, _m{A.at(0).size()}
{
    if(!util::checkRectangular(_A)){
        throw std::invalid_argument("Matrix is not rectangular!");
    }
}

Matrix::Matrix(InitType type, size_t dim) :
_A(Matrix::initMatrix(type, dim)),
_n{_A.size()}, _m{_A.at(0).size()}
{
    if(dim <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }
    if(!util::checkRectangular(_A)){
        throw std::invalid_argument("Matrix is not rectangular!");
    }
}

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

std::unique_ptr<Matrix> Matrix::mult(const util::matrix& B) const {
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
    return std::make_unique<Matrix>(C);
}

std::unique_ptr<Matrix> Matrix::mult(const Matrix& B) const {
    return mult(B.getMatrix());
}

std::unique_ptr<Matrix> Matrix::T() const {
    util::matrix B(_m, util::vec(_n, 0.0));

    for(int i=0; i<_n; i++) {
        for(int j=0; j<_m; j++) {
            B.at(j).at(i) = _A.at(i).at(j);
        }
    }

    return std::make_unique<Matrix>(B);
}

util::vec Matrix::solve(const util::vec& v) const {

}

std::unique_ptr<Matrix> Matrix::inv() const {

}

double Matrix::det() const {

}

util::matrix Matrix::getMatrix() const {
    return _A;
}

void Matrix::print() const {
    util::print(_A);
}

std::pair<double, util::matrix> Matrix::gauss() const {

}