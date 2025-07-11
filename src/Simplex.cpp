#include "Simplex.h"
#include "LinAlg.h"

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

void Simplex::Variable::addConstraint(la::matrix lhs, la::vec rhs, ConstrType type) {
    Simplex::checkSystem(lhs,rhs);
    if(lhs.at(0).size() != _num) {
        throw std::invalid_argument("Matrix doesn't correspond to the variable length!");
    }

    if (_type == unb) {
        // Needs to account for x = x+ - x-
        // Ax = b <=> [A | -A] [x+^t | x-^t]^t = b
        la::conc(lhs, la::scale(lhs, -1.0));
    }

    switch(type) {
        case equal: {
            // Ax = b - canonical form
            // Append A to the end of _lhs, i.e. A_new = [A_old^t | A^t]^t
            la::append(_lhs, std::move(lhs));
            la::append(_rhs, std::move(rhs));
            return;
        };
        case leq: {
            // Add previous slacks
            if(_s_length > 0) {
                la::append(lhs, la::zeros(lhs.size(), _s_length));
            }
            _s_length += lhs.size(); // Number of slack variables for this set of constraints

            // Ax <= b <=> [A | I] [x^t | s^t]^t = b - add slack at the right position
            la::append(lhs, la::eye(lhs.size()));

            la::append(_lhs, std::move(lhs));
            la::append(_rhs, std::move(rhs));
            return;
        };
        case geq: {
            // Ax >= b <=> -Ax <= -b - and proceed as in leq
            la::scale(lhs,-1.0); // Negate the lhs
            la::scale(rhs,-1.0);

            // Add previous slacks
            if(_s_length > 0) {
                la::append(lhs, la::zeros(lhs.size(), _s_length));
            }
            _s_length += lhs.size(); // Number of slack variables for this set of constraints

            // Ax <= b <=> [A | I] [x^t | s^t]^t = b - add slack at the right position
            la::append(lhs, la::eye(lhs.size()));

            la::append(_lhs, std::move(lhs));
            la::append(_rhs, std::move(rhs));
            return;
        };
        default: throw std::invalid_argument("Undefined constraint type!");
    }
}

void Simplex::Variable::addObjective(const la::vec& obj, ObjType type) {
    // Default is min c^t * x
    if(type == min) {
        _obj = obj;
    }
    else {
        _obj = la::scale(obj, -1);
    }

    // x = x+ - x- -> c^t * x = [c^t -c^t] * [x+ x-]
    if(_type == unb) {
        la::append(_obj, la::scale(_obj, -1));
    }
}

std::tuple<la::matrix, la::vec, la::vec> Simplex::Variable::getABC() {
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
    for(la::vec& row : _lhs) {
        if(row.size() < row_length) {
            la::append(row, la::vec(row_length - row.size(), 0.0));
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
            _lhs[i] = la::scale(_lhs[i], -1);
        }
    }

    return {_lhs,_rhs,_obj};
}

la::vec Simplex::Variable::getValue(la::vec sol) {
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

    la::vec result{};
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

void Simplex::addConstraints(const std::string& alias, const la::matrix& lhs, const la::vec& rhs, ConstrType type) {
    std::shared_ptr<Simplex::Variable> var_ptr = findVar(alias);
    if(!var_ptr) {
        throw std::invalid_argument("Alias not assigned to any variable!");
    }
    var_ptr->addConstraint(lhs,rhs,type);
}

void Simplex::addObjective(const std::string& alias, const la::vec& obj, ObjType type) {
    std::shared_ptr<Simplex::Variable> var_ptr = findVar(alias);
    if(!var_ptr) {
        throw std::invalid_argument("Alias not assigned to any variable!");
    }
    var_ptr->addObjective(obj,type);
}

// --- General Solve ---
std::pair<la::vec, Simplex::Status> Simplex::solve() {
    // Formulate the problem
    setup();

    // --- Start Formulate Init
    // Now we have a valid problem formulation
    // Standard form: min c^t*x s.t. Ax=b, x>=0
    int n = _b.size();
    int m = _c.size();

    // Solve the initial feasibility problem >> min s << s.t. [A|I]*[x|s]^t = b, x,s >= 0
    la::matrix A_init = _A;
    la::conc(
        A_init, // Original constraints
        la::eye(n) // Slack variables
    );
    la::vec b_init = _b;
    la::vec c_init = la::vec(m, 0.0);
    la::append(c_init, la::vec(n, 1.0)); // Minimize the slack variables
    la::vec x_init = la::vec(m, 0.0);
    la::vec b_copy = _b;
    la::append(x_init, std::move(b_copy)); // b is the initial solution
    // --- End Formulate Init

    // --- Start Solve Init
    std::pair<la::vec, Simplex::Status> init = solve(std::move(A_init), std::move(b_init), std::move(c_init), std::move(x_init));

    // Initial point : discard slack from the initial solution
    la::vec x = la::subvector(init.first, 0, m);
    la::vec slack = la::subvector(init.first, m, m+n);

    // Initial solution not optimal or slack not 0 <-> No feasible point
    if(init.second != optimal || !std::all_of(slack.begin(), slack.end(), [](const double& val){
        return la::isEqual(val, 0.0);
    })) {
        return {distillSolution(x), infeasible};
    }
    // --- End Solve Init

    std::pair<la::vec, Simplex::Status> final = solve(std::move(_A), std::move(_b), std::move(_c), std::move(x));

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
    if(!la::isRectangular(_A)) {
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

void Simplex::checkSystem(const la::matrix& lhs, const la::vec& rhs) {
    // Must be rectangular, size > 0, matching sizes
    if(!la::isRectangular(lhs)) {
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
        la::append(_b,std::move(b));
        la::append(_c,std::move(c));
        if(iter > 0){
            // Append the 0's
            la::matrix tmp = la::zeros(A.size(), iter);
            la::append(tmp, A); // [0's | A]
            A = std::move(tmp);
        }
        la::append(_A,std::move(A));

        // Current variable starts at iter
        _var_starts.push_back({elem.first, iter});
        // How many new variables added
        iter += c.size();
    });

    // Make _A rectangular - add 0's from the other side
    // iter == _c.size() now, == #variables == #columns of _A
    for(la::vec& row : _A) {
        if(row.size() < iter) {
            la::append(row, la::vec(iter - row.size(), 0.0));
        }
    }

    // Something went wrong
    if(!validSystem()) {
        throw std::runtime_error("Optimization system invalid!");
    }
}

la::vec Simplex::distillSolution(const la::vec& sol_slack) {
    la::vec sol{};

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
        la::append(sol, _vars.at(alias)->getValue(la::subvector(sol_slack, start, end)));
    }
    
    return sol;
}

la::vec Simplex::combineSolution(const la::vec& x_B, const std::vector<int>& B_set, const int& n) {
    la::vec sol(n, 0.0);

    for(int i=0; i<B_set.size(); i++) {
        sol[B_set[i]] = x_B[i];
    }

    return sol;
}

std::pair<la::vec, Simplex::Status> Simplex::solve(la::matrix&& lhs, la::vec&& rhs, la::vec&& obj, la::vec&& init) {
    int maxiter = 1000;

    // -- Problem formulation
    // min c^t*x s.t. Ax=b, x>=0 -> Algorithm for this form

    // --- Start Printing
    std::cout << "--> Initializing Simplex Solver " << std::endl;
    std::cout << "--> Considering problem min(x) c^t*x s.t. A*x=0, x>=0 with: " << std::endl;
    std::cout << "----- A -----" << std::endl;
    la::print(lhs);
    std::cout << "----- b -----" << std::endl;
    la::print(rhs);
    std::cout << "----- c -----" << std::endl;
    la::print(obj);
    std::cout << "----- x -----" << std::endl;
    la::print(init);
    std::cout << "----------" << std::endl;
    // --- End Printing
    
    // ----- Start Initial setting
    std::cout << "--> Initial setting " << std::endl;
    int n = obj.size();
    int m = rhs.size();
    std::cout << "n: " << n << ", m: " << m << std::endl;

    // Get the initial matrices in the *column representation*
    la::matrix lhs_col = la::t(lhs);

    // --- Initialize variables

    // Subsets of the *columns* of the initial LHS matrix
    std::vector<int> N_set; // Normal set -> values of x that are 0
    std::vector<int> B_set; // Basic set -> values of x that are NOT 0 (currently considered)
    N_set.reserve(n-m);
    B_set.reserve(m); // The matrix corresponding to the basic set should be invertible - unique identifier of the corresponding starting edge

    la::matrix a_B{};
    a_B.reserve(m);
    la::matrix a_N{};
    a_N.reserve(n-m);
    la::vec c_N{};
    c_N.reserve(n-m);
    la::vec c_B{};
    c_B.reserve(m);

    // la::vec x_N(n, 0.0); // - always the case
    la::vec x_B{};
    x_B.reserve(m);
    la::vec s_N(n-m, 0.0);

    for(int i=0; i<n; i++){
        if(la::isEqual(init.at(i), 0.0)) {
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
    la::print(N_set);

    std::cout << "----- B_set -----" << std::endl;
    la::print(B_set);

    // --- Initialize Matrices

    // Use matrices with the *column representation*
    la::matrix A_B_t = a_B;
    std::cout << "----- A_B_t -----" << std::endl;
    la::print(A_B_t);
    la::matrix A_N_t = a_N;
    std::cout << "----- A_N_t -----" << std::endl;
    la::print(A_N_t);

    la::matrix A_inv = la::inv(A_B_t);

    int i_N{0}, i_B{0}; // indices to be swapped in each iteration
    // Sets/structures affected:
    // N <-> B
    // c_N <-> c_B
    // A_N_t <-> A_B_t (originally columns, here rows in row representations)
    // s_N and x_B NOT since the counterparts are just 0

    la::vec d{};
    double x_pivot{0.0};

    la::vec lambda{};
    // ----- End Initial setting

    // Progress tracking
    std::vector<double> objective_history;
    const bool show_detailed_output = (n <= 10 && m <= 5); // Only show detailed matrices for small problems
    
    // Calculate initial objective value
    double current_objective = 0.0;
    for(int i = 0; i < B_set.size(); i++) {
        current_objective += c_B[i] * x_B[i];
    }
    objective_history.push_back(current_objective);
    
    std::cout << "--> Starting optimization (initial objective: " << current_objective << ")" << std::endl;
    std::cout << "Progress: ";

    // ----- Start Iterative Steps
    for(int iter=0; iter<maxiter; iter++) {

        // Progress indicator
        if(iter % 10 == 0 || iter < 5) {
            std::cout << "." << std::flush;
        }

        // --- Start Printing (only for small problems or first few iterations)
        if(show_detailed_output || iter < 2) {
            std::cout << "\n-> Iteration " << (iter+1) << " out of " << maxiter 
                      << " (objective: " << current_objective << ")" << std::endl;
            if(show_detailed_output) {
                std::cout << "----- A_B_t -----" << std::endl;
                la::print(A_B_t);
                std::cout << "----- A_N_t -----" << std::endl;
                la::print(A_N_t);
                std::cout << "----- A_inv -----" << std::endl;
                la::print(A_inv);
                std::cout << "----- x_B -----" << std::endl;
                la::print(x_B);
                std::cout << "----- s_N -----" << std::endl;
                la::print(s_N);
                std::cout << "----------" << std::endl;
            }
        }
        // --- End Printing

        // --- Start Iteration calculations
        // 1) l = A_B^(-t)*c_B
        lambda = la::solve(A_B_t, c_B);
        if(show_detailed_output) {
            std::cout << "-- lambda --" << std::endl;
            la::print(lambda);
        }

        // 2) s_N = c_N - A_N^t*l
        s_N = la::sub(c_N, la::mult(A_N_t, lambda));
        if(show_detailed_output) {
            std::cout << "-- s_N --" << std::endl;
            la::print(s_N);
        }
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
            std::cout << "\n--> Optimal solution found!" << std::endl;
            std::cout << "    Final objective value: " << current_objective << std::endl;
            std::cout << "    Iterations completed: " << iter << std::endl;
            
            // Print objective evolution summary
            if(objective_history.size() > 1) {
                std::cout << "    Objective improvement: " << (objective_history[0] - current_objective) << std::endl;
                std::cout << "    Convergence history: ";
                for(size_t i = 0; i < std::min(objective_history.size(), size_t(5)); i++) {
                    std::cout << objective_history[i];
                    if(i < std::min(objective_history.size(), size_t(5)) - 1) std::cout << " -> ";
                }
                if(objective_history.size() > 5) std::cout << " -> ... -> " << current_objective;
                std::cout << std::endl;
            }
            
            return {combineSolution(x_B,B_set,n), optimal};
        }

        // s_N.at(i_N) negative!
        if(show_detailed_output || iter < 2) {
            std::cout << "-> Not optimal yet (entering variable: " << N_set[i_N] << ")" << std::endl;
            std::cout << "i_N: " << i_N << std::endl;
        }
        // --- End Optimality Condition

        // --- Start Pivot
        // i_N: N -> B

        // 4) Let d = A_B^(-1)A_q
        // - If d <= 0 -> problem unbounded. Otherwise consider i : d_i > 0
        d = la::mult(A_inv, la::getRow(A_N_t, i_N)); // A_inv * A_N_t.getRow(i_N) -> d = A_B^(-1) * A_N_t.getRow(i_N)
        // d = A_B_t.T().solve(A_N_t.getRow(i_N));

        if(show_detailed_output) {
            std::cout << "-- d --" << std::endl;
            la::print(d);
        }

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
            std::cout << "\n--> Problem is unbounded!" << std::endl;
            std::cout << "    Direction vector d has no positive components" << std::endl;
            std::cout << "    Iterations completed: " << iter << std::endl;
            return {combineSolution(x_B,B_set,n), unbounded};
        }

        if(show_detailed_output || iter < 2) {
            std::cout << "-> Checked the d condition with" << std::endl;
            std::cout << "x_pivot: " << x_pivot << std::endl;
            std::cout << "i_B: " << i_B << " (leaving variable: " << B_set[i_B] << ")" << std::endl;
        }

        // --- End Pivot

        // --- Start Swapping and Updating
        // 6) Update new x_B = old x_B + Dx_B, swap index p with q in N and B

        // Update x_B
        x_B = la::sub(x_B, la::scale(d, x_pivot));
        x_B.at(i_B) = x_pivot;
        
        // Calculate new objective value
        current_objective = 0.0;
        for(int i = 0; i < B_set.size(); i++) {
            current_objective += c_B[i] * x_B[i];
        }
        objective_history.push_back(current_objective);
        
        if(show_detailed_output || iter < 2) {
            std::cout << "-> Indices to be swapped" << std::endl;
            std::cout << N_set.at(i_N) << ": N -> B " << std::endl;
            std::cout << B_set.at(i_B) << ": B -> N " << std::endl;
            std::cout << "-> New objective value: " << current_objective << std::endl;
        }

        // TODO Sherman Morrison infesible -> try different index from N (matrix not invertible)
        // Sherman-Morrison update for A_inv
        // u = A_N.i_N - A_B.i_B
        // v = [0 ... 1 ... 0], at i_B-th entry
        la::vec v(m, 0.0);
        v.at(i_B) = 1.0;
        la::vec tmp = la::sub(la::getRow(A_N_t, i_N), la::getRow(A_B_t, i_B));
        A_inv = la::ShermanMorrisonInv(A_inv, std::move(tmp), v);

        // Sets/structures affected:
        // N <-> B
        // c_N <-> c_B
        // A_N_t <-> A_B_t (originally columns, here rows in row representations)
        // s_N and x_B NOT since the counterparts are just 0
        std::swap(N_set.at(i_N), B_set.at(i_B));
        std::swap(c_N.at(i_N), c_B.at(i_B));
        la::swapRows(A_N_t, A_B_t, i_N, i_B);
        // --- End Swapping and Updating

    }
    // ----- End Iterative Steps

    std::cout << "\n--> Maximum iterations reached!" << std::endl;
    std::cout << "    Current objective value: " << current_objective << std::endl;
    std::cout << "    Iterations completed: " << maxiter << std::endl;
    
    // Print objective evolution summary
    if(objective_history.size() > 1) {
        std::cout << "    Objective change: " << (objective_history[0] - current_objective) << std::endl;
        std::cout << "    Final convergence: ";
        for(size_t i = std::max(int(objective_history.size()) - 5, 0); i < objective_history.size(); i++) {
            std::cout << objective_history[i];
            if(i < objective_history.size() - 1) std::cout << " -> ";
        }
        std::cout << std::endl;
    }

    throw std::runtime_error("Max. iter count reached and no solution!");
}