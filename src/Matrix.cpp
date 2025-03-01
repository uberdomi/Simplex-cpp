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

Matrix::Matrix(util::matrix A) :
_A(A),
_n{A.size()}, _m{A.at(0).size()}
{}

Matrix::Matrix(InitType type, size_t dim) :
_A(Matrix::initMatrix(type, dim)),
_n{_A.size()}, _m{_A.at(0).size()}
{}