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
        // To be stored in a map of unique_ptrs
        private:
        const int _num; // Once associated with an alias, cannot change its size
        util::matrix _lhs;
        util::vec _rhs;
        util::vec _x_pos, _x_neg, _s; // positive/negative part, slack

        public:
        Variable(int num, VarType var_type);

        void addConstraint(util::matrix lhs, util::vec rhs, ConstrType = equal);

        // [perf] move the own lhs/rhs values and invalidate them -> no longer needed
        // Turn own constraints into the standard form Ax=b
        std::pair<util::matrix, util::vec> getAb();

        // Given the solution vector with x+, x- and slack, return the final variable value
        util::vec getValue(util::vec sol);

    };

    std::unordered_map<std::string, std::unique_ptr<Variable>> _vars;

    public:

    enum Status {
        optimal, infeasible, unbounded
    };

    Simplex()=default;

    // --- New Approach

    void addVariables(const std::string& alias, int num=1, VarType type = pos);

    void addConstraints(const std::string& alias, const util::matrix& lhs, const util::vec& rhs, ConstrType = equal);

    void addObjective(const std::string& alias, const util::vec& obj, ObjType = min);

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
    std::pair<std::pair<util::vec, util::vec>, Simplex::Status> solve();
    // --- Old Approach

    private:
    // --- New Approach
    bool findVar(const std::string& alias);

    // --- Old Approach
    // Helper functions
    bool validSystem() const;

    static std::pair<util::vec, util::vec> distillSolution(const util::vec& x_B, const std::vector<int>& B_set, const int& n);

};

#endif // _SIMPLEX