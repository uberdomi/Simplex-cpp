#ifndef _SIMPLEX
#define _SIMPLEX

#include <list>
#include <map>
#include <memory>
#include <unordered_map>
#include <string>
#include "LinAlg.h"

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

    private:
    // Defined for the final problem formulation
    la::vec _b{}, _c{};
    la::matrix _A{};

    class Variable {
        // [perf] Optimize passing by value/reference - avoid copying data when unnecessary. Generally the outer scope shuld retain its own Ab pairs and Simplex/Variable should use only one copy in total

        // The variable manages its constraints and representation
        // To be stored in a map of shared_ptrs
        private:
        const int _num; // Once associated with an alias, cannot change its size
        const VarType _type;
        la::matrix _lhs;
        la::vec _rhs, _obj;
        int _s_length=0; // Each time an inequality constraint is added, needs to add as many slack variables; the final matrix then needs to be rectangular, with column length _num + _s_length (i.e. add 0's where necessary)

        public:
        Variable(int num, VarType var_type);

        // [pref] Many constraints -> potentially a lot of sparse matrices, inefficient
        void addConstraint(la::matrix lhs, la::vec rhs, ConstrType type = equal);

        void addObjective(const la::vec& obj, ObjType type = min);

        // [perf] move the own lhs/rhs values and invalidate them -> no longer needed
        // Turn own constraints into the standard form Ax=b
        // Make certain b >= 0
        std::tuple<la::matrix, la::vec, la::vec> getABC();

        // Given the solution vector with x+, x- and slack, return the final variable value
        la::vec getValue(la::vec sol);
    };

    // Map for adding variables and in the end aggregating their values
    std::unordered_map<std::string, std::shared_ptr<Variable>> _vars;

    // Stores the start indices for each alias in the final variable vector
    std::vector<std::pair<std::string, int>> _var_starts;

    public:

    enum Status {
        optimal, infeasible, unbounded
    };

    Simplex()=default;

    // --- Optimization

    void addVariables(const std::string& alias, int num=1, VarType type = pos);

    void addConstraints(const std::string& alias, const la::matrix& lhs, const la::vec& rhs, ConstrType type = equal);

    void addObjective(const std::string& alias, const la::vec& obj, ObjType type = min);

    std::pair<la::vec, Simplex::Status> solve();

    private:
    // --- Helper functions

    std::shared_ptr<Variable> findVar(const std::string& alias);

    bool validSystem() const;

    /// @brief Throw an exception if the system of constraints lhs * x = rhs is invalid
    static void checkSystem(const la::matrix& lhs, const la::vec& rhs);

    /// @brief Collect the problem formulations from all stored variables and append them to the final system in the standard form min c^t*x, s.t. Ax=b, x>=0
    void setup();

    /// @brief Solves the system in the standard form min c^t*x, s.t. Ax=b, x>=0, *provided* the valid initial point
    std::pair<la::vec, Simplex::Status> solve(la::matrix&& lhs, la::vec&& rhs, la::vec&& obj, la::vec&& init);

    /// @brief Provided the solution with all reformulations of the variables and slack variables, return the original vector value
    la::vec distillSolution(const la::vec& sol_slack);

    /// @brief The solution with all reformulations (undistilled) is the basic vector filled with 0's at the entries from the normal set
    la::vec combineSolution(const la::vec& x_B, const std::vector<int>& B_set, const int& n);

};

#endif // _SIMPLEX