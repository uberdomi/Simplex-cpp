#ifndef _SIMPLEX
#define _SIMPLEX

#include <list>
#include <map>
#include <memory>
#include <unordered_map>
#include <string>
#include "Utils.h"

enum ConstrType {
    equal, leq, geq
};

enum VarType {
    unb, pos
};

enum ObjType {
    min, max
};

class Simplex {
    // -- Problem formulation
    // min c^t*x s.t. Ax=b, x>=0 -> Algorithm for this form
    // -> For now default formulation Ax<=b and we add slack variables
    private:
    // Defined for the final problem formulation
    util::vec _b{}, _c{};
    util::matrix _A{};

    class Variable {
        // [perf] Optimize passing by value/reference - avoid copying data when unnecessary. Generally the outer scope shuld retain its own Ab pairs and Simplex/Variable should use only one copy in total
        // TODO type, number, constraints (LHS, RHS)
        // The variable manages its constraints and representation
        // To be stored in a map of shared_ptrs
        private:
        const int _num; // Once associated with an alias, cannot change its size
        const VarType _type;
        util::matrix _lhs; // Each row has length _num
        util::vec _rhs, _obj;
        int _s_length=0; // Each time an inequality constraint is added, needs to add as many slack variables; the final matrix then needs to be rectangular, with column length _num + _s_length (i.e. add 0's where necessary)

        public:
        Variable(int num, VarType var_type);

        // [pref] Many constraints -> potentially a lot of sparse matrices, inefficient
        void addConstraint(util::matrix lhs, util::vec rhs, ConstrType constr_type = equal);

        void addObjective(const util::vec& obj, ObjType type = min);

        // [perf] move the own lhs/rhs values and invalidate them -> no longer needed
        // Turn own constraints into the standard form Ax=b
        // Make certain b >= 0
        std::tuple<util::matrix, util::vec, util::vec> getABC();

        // Given the solution vector with x+, x- and slack, return the final variable value
        util::vec getValue(util::vec sol);

    };

    std::unordered_map<std::string, std::shared_ptr<Variable>> _vars;

    // Stores the start indices for each alias in the final variable vector
    std::vector<std::pair<std::string, int>> _var_starts;

    public:

    enum Status {
        optimal, infeasible, unbounded
    };

    Simplex()=default;

    // --- New Approach

    void addVariables(const std::string& alias, int num=1, VarType type = pos);

    void addConstraints(const std::string& alias, const util::matrix& lhs, const util::vec& rhs, ConstrType type = equal);

    void addObjective(const std::string& alias, const util::vec& obj, ObjType type = min);

    std::pair<util::vec, Simplex::Status> solve();

    // --- New Approach

    // --- Old Approach
    // Constraints

    // void addConstraint(const util::vec& lhs, double rhs);
    // void addConstraints(const util::matrix& lhs, const util::vec& rhs);

    // Objective

    // void setMinimization(const util::vec& objective);
    // void setMaximization(const util::vec& objective);

    // Solving

    // {Primal, slack}, status
    // std::pair<std::pair<util::vec, util::vec>, Simplex::Status> solve();
    // --- Old Approach

    private:
    // Helper functions

    // --- Old Approach
    bool validSystem() const;

    // --- New Approach
    std::shared_ptr<Variable> findVar(const std::string& alias);

    /// @brief Throw an exception if the system of constraints lhs * x ? rhs is invalid
    static void checkSystem(const util::matrix& lhs, const util::vec& rhs);

    /// @brief Collect the problem formulations from all stored variables and append them to the final system in the standard form min c^t*x, s.t. Ax=b, x>=0
    void setup();

    /// @brief Solves the system in the standard form min c^t*x, s.t. Ax=b, x>=0, *provided* the valid initial point
    std::pair<util::vec, Simplex::Status> solve(util::matrix&& lhs, util::vec&& rhs, util::vec&& obj, util::vec&& init);

    /// @brief Provided the solution with all reformulations of the variables and slack variables, original vector value
    util::vec distillSolution(const util::vec& sol_slack);

    util::vec distillSolution(const util::vec& x_B, const std::vector<int>& B_set, const int& n);

};

#endif // _SIMPLEX