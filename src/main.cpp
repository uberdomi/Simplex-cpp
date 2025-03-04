#include "Utils.h"
#include "Matrix.h"
#include "Simplex.h"
#include <iostream>

int main(int argc, char** argv) {
    // util::matrix A{{1,2,3}, {2,1,3}, {7,4,2}};
    // util::print(A);

    Matrix A = Matrix({{1,2,3}, {2,1,3}, {-6,-4,2}});
    A.print();

    Matrix B = A.inv();

    std::cout << "--- B = A^(-1) ---" << std::endl;
    B.print();

    std::cout << "--- A * B ---" << std::endl;
    (A * B).print();

    std::cout << "--- B * A ---" << std::endl;
    (B * A).print();

    std::cout << "--- B^(-1)*v = A*v ---" << std::endl;
    util::print(B / util::vec{1,1,1});

    Matrix Ones = Matrix(util::ones(3));

    std::cout << "--- 1_3x3 * 1_3x3 ---" << std::endl;
    Ones.mult(Ones).print();

    std::cout << "--- A + 1_3x3 ---" << std::endl;
    (A + Ones).print();

    std::cout << "--- A | 1_3x3 ---" << std::endl;
    (A | Ones).print();

    std::cout << "--- A + 6*I ---" << std::endl;
    (A + (Matrix(util::eye(3)) * 6)).print();

    std::cout << "--- (A + 6*I)^(-1) ---" << std::endl;
    (A + (Matrix(util::eye(3)) * 6)).inv().print();

    // ----- Simplex Example -----
    std::cout << "----- Simplex Example -----" << std::endl;

    util::matrix A0{{1,3}, {1,1}, {3,1}};
    util::vec b0{15,7,15};
    util::vec c0{1,-1};
    // Considering the LP min(x) c^t*x s.t. A*x<=0, x>=0
    Simplex model{};
    model.addConstraints(A0,b0);
    model.setMaximization(c0);
    auto [sol_slack, status] = model.solve();

    util::print(sol_slack.first);

    util::print(sol_slack.second);

    return 0;
}