#ifndef _SIMPLEX
#define _SIMPLEX

#include <list>
#include "Utils.h"

class Simplex {
    // -- Problem formulation
    // min c^t*x s.t. Ax=b, x>=0 -> Algorithm for this form
    // -> For now default formulation Ax<=b and we add slack variables
    private:
    util::vec _b{}, _c{};
    util::matrix _A{};

    public:

    Simplex()=default;

    enum Status {
        optimal, infeasible, unbounded
    };

    // Constraints

    void addConstraint(const util::vec& coefficients, double rhs);
    void addConstraints(const util::matrix& coefficients, const util::vec& rhs);

    // Objective

    void setMinimization(const util::vec& objective);
    void setMaximization(const util::vec& objective);

    // Solving

    // Primal, slack, status
    std::pair<std::pair<util::vec, util::vec>, Simplex::Status> solve();

    private:
    // Helper functions
    bool validSystem() const;

    static std::pair<util::vec, util::vec> distillSolution(const util::vec& x_B, const std::vector<int>& B_set, const int& n);

};

#endif // _SIMPLEX