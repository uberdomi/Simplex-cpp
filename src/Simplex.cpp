#include "Simplex.h"
#include "Matrix.h"
#include "Utils.h"

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <bits/stdc++.h>

// --- Variable ---

Simplex::Variable::Variable(int num, VarType var_type) : _num(num), _type{var_type}, _lhs{}, _rhs{}, _obj{} {
    if(num < 1) {
        throw std::invalid_argument("Variable length must be greater than 0!");
    }
}

void Simplex::Variable::addConstraint(util::matrix lhs, util::vec rhs, ConstrType type) {
    Simplex::checkSystem(lhs,rhs);
    if(lhs.at(0).size() != _num) {
        throw std::invalid_argument("Matrix doesn't correspond to the variable length!");
    }

    if (_type == unb) {
        // Needs to account for x = x+ - x-
        // Ax = b <=> [A | -A] [x+^t | x-^t]^t = b
        lhs = (Matrix(lhs) | (Matrix(lhs) * (-1))).get();
    }

    switch(type) {
        case equal: {
            // Ax = b - canonical form
            // Append A to the end of _lhs, i.e. A_new = [A_old^t | A^t]^t
            util::append(_lhs, std::move(lhs));
            util::append(_rhs, std::move(rhs));
            return;
        };
        case leq: {
            // Add previous slacks
            if(_s_length > 0) {
                lhs = (Matrix(lhs) | Matrix(util::zeros(lhs.size(), _s_length))).get();
            }
            _s_length += lhs.size(); // Number of slack variables for this set of constraints

            // Ax <= b <=> [A | I] [x^t | s^t]^t = b - add slack at the right position
            lhs = (Matrix(lhs) | Matrix(util::eye(lhs.size()))).get();

            util::append(_lhs, std::move(lhs));
            util::append(_rhs, std::move(rhs));
            return;
        };
        case geq: {
            // Ax >= b <=> -Ax <= -b - and proceed as in leq
            lhs = (Matrix(lhs) * (-1)).get();
            rhs = util::scale(rhs, -1);

            // Add previous slacks
            if(_s_length > 0) {
                lhs = (Matrix(lhs) | Matrix(util::zeros(lhs.size(), _s_length))).get();
            }
            _s_length += lhs.size(); // Number of slack variables for this set of constraints

            // Ax <= b <=> [A | I] [x^t | s^t]^t = b - add slack at the right position
            lhs = (Matrix(lhs) | Matrix(util::eye(lhs.size()))).get();

            util::append(_lhs, std::move(lhs));
            util::append(_rhs, std::move(rhs));
            return;
        };
        default: throw std::invalid_argument("Undefined constraint type!");
    }
}

void Simplex::Variable::addObjective(const util::vec& obj, ObjType type) {
    // Default is min c^t * x
    if(type == min) {
        _obj = obj;
    }
    else {
        _obj = util::scale(obj, -1);
    }

    // x = x+ - x- -> c^t * x = [c^t -c^t] * [x+ x-]
    if(_type == unb) {
        util::append(_obj, util::scale(_obj, -1));
    }
}

std::tuple<util::matrix, util::vec, util::vec> Simplex::Variable::getABC() {
    // Add 0's where slack is present, to bring A to a rectangular form
    // Length corresponds to the variable vector [x+ x- s]
    int x_length{};
    if(_type == unb) {
        x_length = 2*_num;
    }
    else {
        x_length = _num;
    }
    int row_length = x_length + _s_length;
    for(util::vec& row : _lhs) {
        if(row.size() < row_length) {
            util::append(row, util::vec(row_length - row.size(), 0.0));
        }
    }

    // Add slack to the objective
    _obj.reserve(row_length);
    for(int i=_obj.size(); i<_obj.capacity(); i++) {
        _obj[i] = 0.0;
    }

    // Make sure rhs is positive
    for(int i=0; i<_rhs.size(); i++) {
        if(_rhs[i] < 0) {
            _rhs[i] = -_rhs[i];
            _lhs[i] = util::scale(_lhs[i], -1);
        }
    }

    return {_lhs,_rhs,_obj};
}

util::vec Simplex::Variable::getValue(util::vec sol) {
    // Solution has length of row_length, is of form [x+ x- s]
    int x_length{};
    if(_type == unb) {
        x_length = 2*_num;
    }
    else {
        x_length = _num;
    }
    int row_length = x_length + _s_length;

    if(sol.size() != row_length) {
        throw std::invalid_argument("Solution doesn't correspond to the internal variable parameters");
    }

    util::vec result{};
    result.reserve(_num);
    if(_type == unb) {
        // x = x+ - x-
        for(int i=0; i<_num; i++) {
            result[i] = sol.at(i) - sol.at(_num+i);
        }
    }
    else {
        // x = x+
        std::transform(sol.begin(), sol.begin()+_num, result.begin(),
        [](double val){return val;});
    }

    return result;
}

// --- Optimization ---

void Simplex::addVariables(const std::string& alias, int num, VarType type){
    if(findVar(alias)) {
        throw std::invalid_argument("Alias already used for the variable!");
    }

    _vars.insert({alias, std::make_shared<Simplex::Variable>(num, type)});
}

void Simplex::addConstraints(const std::string& alias, const util::matrix& lhs, const util::vec& rhs, ConstrType type) {
    std::shared_ptr<Simplex::Variable> var_ptr = findVar(alias);
    if(!var_ptr) {
        throw std::invalid_argument("Alias not assigned to any variable!");
    }
    var_ptr->addConstraint(lhs,rhs,type);
}

void Simplex::addObjective(const std::string& alias, const util::vec& obj, ObjType type) {
    std::shared_ptr<Simplex::Variable> var_ptr = findVar(alias);
    if(!var_ptr) {
        throw std::invalid_argument("Alias not assigned to any variable!");
    }
    var_ptr->addObjective(obj,type);
}

// --- General Solve ---
std::pair<util::vec, Simplex::Status> Simplex::solve() {
    // Formulate the problem
    setup();

    // --- Start Formulate Init
    // Now we have a valid problem formulation
    // Standard form: min c^t*x s.t. Ax=b, x>=0
    int n = _b.size();
    int m = _c.size();

    // Solve the initial feasibility problem >> min s << s.t. [A|I]*[x|s]^t = b, x,s >= 0
    util::matrix A_init = (Matrix(_A) | util::eye(n)).get();
    util::vec b_init = _b;
    util::vec c_init = util::vec(m, 0.0);
    util::append(c_init, util::vec(n, 1.0)); // Minimize the slack variables
    util::vec x_init = util::vec(m, 0.0);
    util::vec b_copy = _b;
    util::append(x_init, std::move(b_copy)); // b is the initial solution
    // --- End Formulate Init

    // --- Start Solve Init
    std::pair<util::vec, Simplex::Status> init = solve(std::move(A_init), std::move(b_init), std::move(c_init), std::move(x_init));

    // Initial point : discard slack from the initial solution
    util::vec x = util::subvector(init.first, 0, m);
    util::vec slack = util::subvector(init.first, m, m+n);

    // Initial solution not optimal or slack not 0 <-> No feasible point
    if(init.second != optimal || !std::all_of(slack.begin(), slack.end(), [](const double& val){
        return util::isEqual(val, 0.0);
    })) {
        return {distillSolution(x), infeasible};
    }
    // --- End Solve Init

    std::pair<util::vec, Simplex::Status> final = solve(std::move(_A), std::move(_b), std::move(_c), std::move(x));

    return {distillSolution(final.first), final.second};
}

// --- Helper functions ---

std::shared_ptr<Simplex::Variable> Simplex::findVar(const std::string& alias) {
    auto it = _vars.find(alias);
    if(_vars.find(alias) == _vars.end()) {
        return nullptr;
    }
    else {
        return _vars.at(alias);
    }
}

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

    // Check whether b >= 0
    return std::all_of(_b.begin(), _b.end(), [](const double& val){
        return val >= 0;
    });
}

void Simplex::checkSystem(const util::matrix& lhs, const util::vec& rhs) {
    // Must be rectangular, size > 0, matching sizes
    if(!util::isRectangular(lhs)) {
        throw std::invalid_argument("Matrix not rectangular!");
    }
    if(lhs.size() == 0){
        throw std::invalid_argument("Empty constraint!");
    }
    if(lhs.size() != rhs.size()){
        throw std::invalid_argument("Constraint dimensions aren't the same!");
    }
}

void Simplex::setup() {
    // Concatenate all the problem formulations into the final system
    int iter=0;
    std::for_each(_vars.begin(), _vars.end(), [this, &iter](const std::pair<std::string, std::shared_ptr<Variable>>& elem) {
        auto [A,b,c] = elem.second->getABC();
        util::append(_b,std::move(b));
        util::append(_c,std::move(c));
        if(iter > 0){
            // Append the 0's
            A = (Matrix(util::zeros(A.size(), iter)) | A).get();
        }
        util::append(_A,std::move(A));

        // Current variable starts at iter
        _var_starts.push_back({elem.first, iter});
        // How many new variables added
        iter += c.size();
    });

    // Make _A rectangular - add 0's from the other side
    // iter == _c.size() now, == #variables == #columns of _A
    for(util::vec& row : _A) {
        if(row.size() < iter) {
            util::append(row, util::vec(iter - row.size(), 0.0));
        }
    }

    // Something went wrong
    if(!validSystem()) {
        throw std::runtime_error("Optimization system invalid!");
    }
}

util::vec Simplex::distillSolution(const util::vec& sol_slack) {
    util::vec sol{};

    int start=0;
    int end=0;
    std::string alias{};
    for (int i=0; i<_var_starts.size(); i++) {
        alias = _var_starts.at(i).first;
        start = _var_starts.at(i).second;
        if(i+1 < _var_starts.size()) {
            end = start = _var_starts.at(i+1).second;
        }
        else {
            end = sol_slack.size();
        }

        // Distill the value from the subvector of the solution vecotr, corresponding to the variable associated with the alias
        util::append(sol, _vars.at(alias)->getValue(util::subvector(sol_slack, start, end)));
    }
    
    return sol;
}

util::vec Simplex::combineSolution(const util::vec& x_B, const std::vector<int>& B_set, const int& n) {
    util::vec sol(n, 0.0);

    for(int i=0; i<B_set.size(); i++) {
        sol[B_set[i]] = x_B[i];
    }

    return sol;
}

std::pair<util::vec, Simplex::Status> Simplex::solve(util::matrix&& lhs, util::vec&& rhs, util::vec&& obj, util::vec&& init) {
    int maxiter = 1000;

    // -- Problem formulation
    // min c^t*x s.t. Ax=b, x>=0 -> Algorithm for this form

    // --- Start Printing
    std::cout << "--> Initializing Simplex Solver " << std::endl;
    std::cout << "--> Considering problem min(x) c^t*x s.t. A*x=0, x>=0 with: " << std::endl;
    std::cout << "----- A -----" << std::endl;
    util::print(lhs);
    std::cout << "----- b -----" << std::endl;
    util::print(rhs);
    std::cout << "----- c -----" << std::endl;
    util::print(obj);
    std::cout << "----- x -----" << std::endl;
    util::print(init);
    std::cout << "----------" << std::endl;
    // --- End Printing
    
    // ----- Start Initial setting
    std::cout << "--> Initial setting " << std::endl;
    int n = obj.size();
    int m = rhs.size();
    std::cout << "n: " << n << ", m: " << m << std::endl;

    // Get the initial matrices in the *column representation*
    util::matrix lhs_col = Matrix(lhs).T().get(); // now each row corresponds to the column

    // --- Initialize variables

    // Subsets of the *columns* of the initial LHS matrix
    std::vector<int> N_set; // Normal set -> values of x that are 0
    std::vector<int> B_set; // Basic set -> values of x that are NOT 0 (currently considered)
    N_set.reserve(n-m);
    B_set.reserve(m); // The matrix corresponding to the basic set should be invertible - unique identifier of the corresponding starting edge

    util::matrix a_B{};
    a_B.reserve(m);
    util::matrix a_N{};
    a_N.reserve(n-m);
    util::vec c_N{};
    c_N.reserve(n-m);
    util::vec c_B{};
    c_B.reserve(m);

    // util::vec x_N(n, 0.0); // - always the case
    util::vec x_B{};
    x_B.reserve(m);
    util::vec s_N(n-m, 0.0);

    for(int i=0; i<n; i++){
        if(util::isEqual(init.at(i), 0.0)) {
            N_set.push_back(i);
            a_N.push_back(lhs_col.at(i));
            c_N.push_back(obj.at(i));
        }
        else {
            B_set.push_back(i);
            a_B.push_back(lhs_col.at(i));
            c_B.push_back(obj.at(i));
            x_B.push_back(init.at(i));
        }
    }

    // TODO : more than m indices can be selected -> needs to select a valid subset
    if(B_set.size() != m) {
        throw std::runtime_error("Cannot select a valid initial basic set!");
    }

    std::cout << "----- N_set -----" << std::endl;
    util::print(N_set);

    std::cout << "----- B_set -----" << std::endl;
    util::print(B_set);

    // --- Initialize Matrices

    // Use matrices with the *column representation*
    Matrix A_B_t = Matrix(a_B);
    std::cout << "----- A_B_t -----" << std::endl;
    A_B_t.print();
    Matrix A_N_t = Matrix(a_N);
    std::cout << "----- A_N_t -----" << std::endl;
    A_N_t.print();

    Matrix A_inv = A_B_t.inv();

    int i_N{0}, i_B{0}; // indices to be swapped in each iteration
    // Sets/structures affected:
    // N <-> B
    // c_N <-> c_B
    // A_N_t <-> A_B_t (originally columns, here rows in row representations)
    // s_N and x_B NOT since the counterparts are just 0

    util::vec d{};
    double x_pivot{0.0};

    util::vec lambda{};
    // ----- End Initial setting

    // ----- Start Iterative Steps
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
        // 1) l = A_B^(-t)*c_B
        lambda = A_B_t.solve(c_B);
        std::cout << "-- lambda --" << std::endl;
        util::print(lambda);

        // 2) s_N = c_N - A_N^t*l
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
            return {combineSolution(x_B,B_set,n), optimal};
        }

        // s_N.at(i_N) negative!
        std::cout << "-> Not optimal yet" << std::endl;
        std::cout << "i_N: " << i_N << std::endl;
        // --- End Optimality Condition

        // --- Start Pivot
        // i_N: N -> B

        // 4) Let d = A_B^(-1)A_q
        // - If d <= 0 -> problem unbounded. Otherwise consider i : d_i > 0
        d = A_inv * A_N_t.getRow(i_N);
        // d = A_B_t.T().solve(A_N_t.getRow(i_N));

        std::cout << "-- d --" << std::endl;
        util::print(d);

        // 5) Select p = argmin i {x_i / d_i} - where d_i > 0 and i from B (to be swapped away)
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
            return {combineSolution(x_B,B_set,n), unbounded};
        }

        std::cout << "-> Checked the d condition with" << std::endl;
        std::cout << "x_pivot: " << x_pivot << std::endl;
        std::cout << "i_B: " << i_B << std::endl;

        // --- End Pivot

        // --- Start Swapping and Updating
        // 6) Update new x_B = old x_B + Dx_B, swap index p with q in N and B

        // Update x_B
        x_B = util::sub(x_B, util::scale(d, x_pivot));
        x_B.at(i_B) = x_pivot;
        
        std::cout << "-> Indices to be swapped" << std::endl;
        std::cout << N_set.at(i_N) << ": N -> B " << std::endl;
        std::cout << B_set.at(i_B) << ": B -> N " << std::endl;

        // TODO Sherman Morrison infesible -> try different index from N (matrix not invertible)
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
    // ----- End Iterative Steps

    throw std::runtime_error("Max. iter count reached and no solution!");
}