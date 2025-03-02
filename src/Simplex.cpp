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
    for(int i=0; i<objective.size(); i++) {
        _c.push_back(-objective.at(i));
    }
}

// Solving

std::pair<util::vec, Simplex::Status> Simplex::solve(){
    if(!validSystem()) {
        throw std::runtime_error("Optimization system invalid!");
    }

    int maxiter = 1000;

    // -- Problem formulation
    // min c^t*x s.t. Ax=b, x>=0 -> Algorithm for this form
    // -> For now default formulation Ax<=b and we add slack variables

    // 0) TODO [later] Find feasible point
    // for now assume [0, ..., 0, b^t] is the initial *feasible* point

    // --- Start Printing
    std::cout << "--> Initializing Simplex Solver " << std::endl;
    std::cout << "--> Considering problem min(x) c^t*x s.t. A*x<=0, x>=0 with: " << std::endl;
    std::cout << "----- A -----" << std::endl;
    util::print(_A);
    std::cout << "----- b -----" << std::endl;
    util::print(_b);
    std::cout << "----- c -----" << std::endl;
    util::print(_c);
    std::cout << "----------" << std::endl;
    // --- End Printing
    
    // --- Start Initial setting
    std::cout << "--> Initial setting " << std::endl;
    int m = _b.size();
    int n = _c.size();
    std::cout << "n: " << n << ", m: " << m << std::endl;

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
    std::cout << "----- A_B_t -----" << std::endl;
    A_B_t.print();
    Matrix A_N_t = Matrix(_A).T();
    std::cout << "----- A_N_t -----" << std::endl;
    A_N_t.print();

    util::vec c_N = _c;
    util::vec c_B(m, 0.0);

    // util::vec x_N(n, 0.0); // - always the case
    util::vec x_B = _b;
    util::vec s_N{};

    int p{0}, p_i{0}, q{0}; // indices to be swapped in each iteration
    util::vec d{};
    double x_q{0.0};

    util::vec lambda{};
    // --- End Initial setting

    for(int iter=0; iter<maxiter; iter++) {

        // --- Start Printing
        std::cout << "-> Iteration " << (iter+1) << " out of " << maxiter << std::endl;
        std::cout << "----- A_B_t -----" << std::endl;
        A_B_t.print();
        std::cout << "----- A_N_t -----" << std::endl;
        A_N_t.print();
        std::cout << "----- x_B -----" << std::endl;
        util::print(x_B);
        std::cout << "----- s_N -----" << std::endl;
        util::print(s_N);
        std::cout << "----------" << std::endl;
        // --- End Printing

        // Iteration calculations
        lambda = A_B_t.solve(c_B);
        std::cout << "-- lambda --" << std::endl;
        util::print(lambda);

        s_N = util::sub(c_N, A_N_t * lambda);
        std::cout << "-- s_N --" << std::endl;
        util::print(s_N);

        // TODO check s_N optimality
        if(std::all_of(s_N.begin(), s_N.end(), [](const double& val) {
            return val >= 0;
        })) {
            // TODO return x corresponding to B and N
            std::cout << "--> Exiting: problem optimal!" << std::endl;
            return {discardSlack(x_B,B_set), optimal};
        }

        std::cout << "-> Not optimal yet" << std::endl;

        // --- Start Pivot
        // Pick indices to be swapped away
        q = N_set.front();
        std::cout << "-- q --" << std::endl;
        std::cout << q << std::endl;

        d = A_B_t.T().solve(A_N_t.getRow(0));
        std::cout << "-- d --" << std::endl;
        util::print(d);

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
        if(util::isEqual(x_q, -1.0)) {
            // TODO check d condition
            // <=> d <= 0 -> unbounded
            std::cout << "--> Exiting: problem unbounded!" << std::endl;
            return {discardSlack(x_B,B_set), unbounded};
        }

        std::cout << "-> Checked the d condition with" << std::endl;
        std::cout << "x_q: " << x_q << std::endl;
        std::cout << "p_i: " << p_i << std::endl;

        // // Get the p_i'th element of the list
        // it = B_set.begin();
        // std::advance(it, p_i);
        // p = *it;
        // Or with that
        it = std::next(B_set.begin(), p_i);  // Move iterator to the p_i-th element

        // Update elements
        std::cout << "-> Updating elements" << std::endl;
        util::print(x_B);
        util::print(util::scale(d, x_q));

        x_B = util::add(x_B, util::scale(d, x_q));

        std::cout << "-> Swapping rows" << std::endl;
        std::cout << "p: " << p << std::endl;
        std::cout << "q: " << q << std::endl;

        // TODO swap indices q, p from N, B resp.
        // From A_N we always pop from the top
        // In A_B we explicitly found the index p_i to replace the row with the row from A_N
        Matrix::swapRows(A_N_t, A_B_t, 0, p_i);
        // Newly swapped-out row emplaced back
        A_N_t.swapToBack(0);

        // Also swap c
        std::swap(c_B.at(p_i), c_N.at(0));
        c_N.push_back(std::move(c_N.front()));  // Move first row to the back
        c_N.erase(c_N.begin());  // Remove the now-empty first row

        // N: +p -q
        // B: -p +q
        N_set.pop();
        N_set.push(p);
        *it = q;

        // --- End Pivot

    }

    throw std::runtime_error("Max. iter count reached and no solution!");
}

// Helper functions

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

util::vec Simplex::discardSlack(const util::vec& x_B, const std::list<int>& B_set){
    int n = x_B.size();
    util::vec x(n, 0.0);
    for(const int& i : B_set) {
        if(i < n) {
            x.at(i) = x_B.at(i);
        }
    }

    return x;
}