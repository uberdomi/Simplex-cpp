#include "Simplex.h"
#include "Matrix.h"
#include "Utils.h"

#include <stdexcept>
#include <algorithm>
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

std::pair<std::pair<util::vec, util::vec>, Simplex::Status> Simplex::solve(){
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

    std::vector<int> N_set(n); // Normal set -> starting with regular variables
    std::vector<int> B_set(m); // Basic set -> starting with slack variables

    // Assign increasing values starting from 0
    std::iota(N_set.begin(), N_set.end(), 0);

    // Assign increasing values starting from n
    std::iota(B_set.begin(), B_set.end(), n);

    std::cout << "----- N_set -----" << std::endl;
    util::print(N_set);

    std::cout << "----- B_set -----" << std::endl;
    util::print(B_set);

    // Use matrices with the *column representation*
    Matrix A_B_t = Matrix(util::eye(m));
    std::cout << "----- A_B_t -----" << std::endl;
    A_B_t.print();
    Matrix A_N_t = Matrix(_A).T();
    std::cout << "----- A_N_t -----" << std::endl;
    A_N_t.print();

    Matrix A_inv = A_B_t.inv();

    util::vec c_N = _c;
    util::vec c_B(m, 0.0);

    // util::vec x_N(n, 0.0); // - always the case
    util::vec x_B = _b;
    util::vec s_N(n, 0.0);

    int i_N{0}, i_B{0}; // indices to be swapped in each iteration
    // Sets/structures affected:
    // N <-> B
    // c_N <-> c_B
    // A_N_t <-> A_B_t (originally columns, here rows in row representations)
    // s_N and x_B NOT since the counterparts are just 0

    util::vec d{};
    double x_pivot{0.0};

    util::vec lambda{};
    // --- End Initial setting

    for(int iter=0; iter<maxiter; iter++) {

        // --- Start Printing
        std::cout << "-> Iteration " << (iter+1) << " out of " << maxiter << std::endl;
        std::cout << "----- A_B_t -----" << std::endl;
        A_B_t.print();
        std::cout << "----- A_N_t -----" << std::endl;
        A_N_t.print();
        std::cout << "----- A_inv -----" << std::endl;
        A_inv.print();
        std::cout << "----- x_B -----" << std::endl;
        util::print(x_B);
        std::cout << "----- s_N -----" << std::endl;
        util::print(s_N);
        std::cout << "----------" << std::endl;
        // --- End Printing

        // --- Start Iteration calculations
        lambda = A_B_t.solve(c_B);
        std::cout << "-- lambda --" << std::endl;
        util::print(lambda);

        s_N = util::sub(c_N, A_N_t * lambda);
        std::cout << "-- s_N --" << std::endl;
        util::print(s_N);
        // --- End Iteration calculations

        // --- Start Optimality Condition
        // Find the first negative entry of s_N
        i_N = -1;
        for(int i=0; i<n; i++) {
            if(s_N.at(i) < 0){
                i_N = i;
                break;
            }
        }
        if(i_N == -1) {
            // s_N non-negative -> optimal solution found
            std::cout << "--> Exiting: problem optimal!" << std::endl;
            return {distillSolution(x_B,B_set,n), optimal};
        }

        // s_N.at(i_N) negative!
        std::cout << "-> Not optimal yet" << std::endl;
        std::cout << "i_N: " << i_N << std::endl;

        // --- End Optimality Condition

        // --- Start Pivot
        // i_N: N -> B

        d = A_inv * A_N_t.getRow(i_N);
        // d = A_B_t.T().solve(A_N_t.getRow(i_N));

        std::cout << "-- d --" << std::endl;
        util::print(d);

        // TODO check d condition
        // TODO iterate over positive entries of d
        i_B = 0;
        x_pivot = -1.0; // All entries should be positive -> x_pivot < 0 <=> not chosen yet
        for(int i=0; i<m; i++) {
            if(d.at(i) > 0.0){
                if(x_pivot < 0.0) {
                    i_B = i;
                    x_pivot = x_B.at(i) / d.at(i);
                }
                if(x_B.at(i) / d.at(i) < x_pivot) {
                    i_B = i;
                    x_pivot = x_B.at(i) / d.at(i);
                }
            }
        }
        if(x_pivot < 0) {
            // <=> d <= 0 -> unbounded
            std::cout << "--> Exiting: problem unbounded!" << std::endl;
            return {distillSolution(x_B,B_set,n), unbounded};
        }

        std::cout << "-> Checked the d condition with" << std::endl;
        std::cout << "x_pivot: " << x_pivot << std::endl;
        std::cout << "i_B: " << i_B << std::endl;

        // --- End Pivot

        // --- Start Swapping and Updating

        // Update x_B
        x_B = util::sub(x_B, util::scale(d, x_pivot));
        x_B.at(i_B) = x_pivot;
        
        std::cout << "-> Indices to be swapped" << std::endl;
        std::cout << N_set.at(i_N) << ": N -> B " << std::endl;
        std::cout << B_set.at(i_B) << ": B -> N " << std::endl;

        // Sherman-Morrison update for A_inv
        // u = A_N.i_N - A_B.i_B
        // v = [0 ... 1 ... 0], at i_B-th entry
        util::vec v(m, 0.0);
        v.at(i_B) = 1.0;
        A_inv = std::move(A_inv.ShermanMorrison(util::sub(A_N_t.getRow(i_N), A_B_t.getRow(i_B)), v));

        // Sets/structures affected:
        // N <-> B
        // c_N <-> c_B
        // A_N_t <-> A_B_t (originally columns, here rows in row representations)
        // s_N and x_B NOT since the counterparts are just 0
        std::swap(N_set.at(i_N), B_set.at(i_B));
        std::swap(c_N.at(i_N), c_B.at(i_B));
        Matrix::swapRows(A_N_t, A_B_t, i_N, i_B);

        // --- End Swapping and Updating

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

std::pair<util::vec, util::vec> Simplex::distillSolution(const util::vec& x_B, const std::vector<int>& B_set, const int& n){
    int m = B_set.size();

    util::vec sol(n, 0.0);
    util::vec slack(m, 0.0);

    // i-th entry of x_B corresponds to the real B_set(i)-th entry
    for(int i=0; i<m; i++) {
        if(B_set.at(i) >= n) {
            slack.at(B_set.at(i) - n) = x_B.at(i);
        }
        else {
            sol.at(B_set.at(i)) = x_B.at(i);
        }
    }

    return {sol,slack};
}