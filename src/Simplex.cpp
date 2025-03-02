#include "Simplex.h"
#include "Matrix.h"
#include "Utils.h"

#include <stdexcept>
#include <algorithm>
#include <queue>
#include <set>
#include <list>
#include <iostream>
#include <bits/stdc++.h> 

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

    // -- Problem formulation
    // min c^t*x s.t. Ax=b, x>=0 -> Algorithm for this form
    // -> For now default formulation Ax<=b and we add slack variables

    // 0) [later] Find feasible point
    // for now assume [0, ..., 0, b^t] is the initial *feasible* point
    
    // --- Initial setting ---
    int m = _b.size();
    int n = _c.size();

    std::queue<int> N_set{}; // Normal set
    std::list<int> B_set{}; // Basic set
    std::list<int>::iterator it = B_set.begin();

    for(int i=0; i<n; i++) {
        N_set.push(i);
    }

    for(int i=n; i<m+n; i++) {
        B_set.push_back(i);
    }

    // Use matrices with the *column representation*
    Matrix A_B_t = Matrix(util::eye(m));
    Matrix A_N_t = Matrix(_A).T();
    Matrix A_full_t = A_N_t | A_B_t;

    util::vec c_N = _c;
    util::vec c_B(m, 0.0);

    // util::vec x_N(n, 0.0); // - always the case
    util::vec x_B = _b;
    util::vec s_N{};

    int p{0}, p_i{0}, q{0}; // indices to be swapped in each iteration
    util::vec d{};
    double x_q{0.0};

    util::vec lambda{};
    // Iteration calculations
    lambda = A_B_t.solve(c_B);
    s_N = util::sub(c_N, A_N_t * lambda);
    // TODO check s_N optimality
    // TODO return x corresponding to B and N

    // --- Start Pivot
    // Pick indices to be swapped away
    q = N_set.front();
    d = A_B_t.T().solve(A_full_t.getRow(q));

    // TODO check d condition
    // TODO iterate over positive entries of d
    p_i = 0;
    x_q = -1.0; // All entries should be positive -> x_q < 0 <=> not chosen yet
    for(int i=0; i<m; i++) {
        if(d.at(i) > 0.0){
            if(x_q < 0.0) {
                p_i = i;
                x_q = x_B.at(i) / d.at(i);
            }
            if(x_B.at(i) / d.at(i) < x_q) {
                p_i = i;
                x_q = x_B.at(i) / d.at(i);
            }
        }
    }
    // Get the p_i'th element of the list
    it = B_set.begin();
    std::advance(it, p_i);
    p = *it;

    // Update elements
    x_B = util::add(x_B, util::scale(d, x_q));

    // TODO swap indices q, p from N, B resp.

    // --- End Pivot

    // Printing steps


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