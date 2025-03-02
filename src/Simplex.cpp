#include "Simplex.h"
#include "Matrix.h"
#include "Utils.h"

#include <stdexcept>
#include <algorithm>

// Constraints

void Simplex::addConstraint(const util::vec& coefficients, double rhs){
    // a_i^t * x <= b_i 
    _A.push_back(coefficients);
    _b.push_back(rhs);
}

void Simplex::addConstraints(const util::matrix& coefficients, const util::vec& rhs){
    // A * x <= b

    if(coefficients.size() != rhs.size()) {
        throw std::invalid_argument("Incompatible sizes of coefficient matrix and the right hand side!");
    }

    for(int i=0; i<rhs.size(); i++) {
        addConstraint(coefficients.at(i), rhs.at(i));
    }
}

// Objective

void Simplex::setMinimization(const util::vec& objective){
    // min c*x
    _c = objective;
}

void Simplex::setMaximization(const util::vec& objective){
    // max c*x <=> min -c*x
    _c.reserve(objective.size());
    std::transform(objective.begin(), objective.end(), _c.begin(), 
    [](const double& val) {return -val;}
    );
}

// Solving

std::pair<util::vec, Simplex::Status> Simplex::solve(){
    if(!validSystem()) {
        throw std::runtime_error("Optimization system invalid!");
    }
}

// System validation

bool Simplex::validSystem() const {
    if(!util::isRectangular(_A)) {
        return false;
    }
    if(_A.size() != _b.size()) {
        return false;
    }

    if(_A.size() > 0 && _A.at(0).size() != _c.size()) {
        return false;
    }

    return true;
}