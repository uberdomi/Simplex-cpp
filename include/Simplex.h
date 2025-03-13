#ifndef _SIMPLEX
#define _SIMPLEX

#include <list>
#include <map>
#include <unordered_map>
#include <string>
#include "Utils.h"

class Simplex {
    // -- Problem formulation
    // min c^t*x s.t. Ax=b, x>=0 -> Algorithm for this form
    // -> For now default formulation Ax<=b and we add slack variables
    private:
    // Defined for the final problem formulation
    util::vec _b{}, _c{};
    util::matrix _A{};

    class Variable {
        // TODO type, number, constraints (LHS, RHS)
        // To be stored in a map of unique_ptrs
    };

    std::unordered_map<std::string, Variable> _vars;

    public:

    enum ConstrType {
        equal, leq, geq
    };

    enum VarType {
        unb, pos
    };

    Simplex()=default;

    enum Status {
        optimal, infeasible, unbounded
    };

    // --- New Approach

    void AddConstraints(const int number, const std::string& alias, const VarType& type = pos);

    // --- New Approach

    // --- Old Approach
    // Constraints

    void addConstraint(const util::vec& coefficients, double rhs);
    void addConstraints(const util::matrix& coefficients, const util::vec& rhs);

    // Objective

    void setMinimization(const util::vec& objective);
    void setMaximization(const util::vec& objective);

    // Solving

    // {Primal, slack}, status
    std::pair<std::pair<util::vec, util::vec>, Simplex::Status> solve();
    // --- Old Approach

    private:
    // Helper functions
    bool validSystem() const;

    static std::pair<util::vec, util::vec> distillSolution(const util::vec& x_B, const std::vector<int>& B_set, const int& n);

};

#endif // _SIMPLEX