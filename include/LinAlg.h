#ifndef _LINALG
#define _LINALG

#include <vector>

namespace la {
    // Helper structures
    using vec = std::vector<double>;
    using matrix = std::vector<vec>;

    // Comparisons
    bool isEqual(const double& x, const double& y, const double& tol=1e-6);
    bool isEqual(const vec& v, const vec& w, const double& tol=1e-6);
    bool isEqual(const matrix& A, const matrix& B, const double& tol=1e-6);

    bool isRectangular(const matrix& A);

    // Vector manipulations
    void scaleVec(vec& v, const double& sf);
    void subtractVec(vec& v, const vec& w, const double& sf=1.0);

    // Printing
    void print(const vec& v);
    void print(const std::vector<int>& v);
    void print(const matrix& A);

    // Conversions
    matrix vecAsMatrix(const vec& v);
    vec matrixAsVec(const matrix& A);

    // Common matrices
    matrix eye(size_t n);
    matrix ones(size_t n, size_t m=0);
    matrix zeros(size_t n, size_t m=0);

    matrix t(const matrix& A);

    // Glueing matrices
    void append(matrix& A, const matrix& B); // [A^t | B^t]^t
    matrix operator^(const matrix& A, const matrix& B);
    matrix& operator^=(matrix& A, const matrix& B);

    void append(matrix& A, const vec& v); // Adding new row
    matrix operator^(const matrix& A, const vec& v);
    matrix& operator^=(matrix& A, const vec& v);

    void conc(vec& v, const vec& w);
    vec operator|(const vec& v, const vec& w);
    vec& operator|=(vec& v, const vec& w);

    void conc(matrix& A, const matrix& B); // [A | B]
    matrix operator|(const matrix& A, const matrix& B);
    matrix& operator|=(matrix& A, const matrix& B);

    void conc(matrix& A, const vec& v); // Adding to each row
    matrix operator|(const matrix& A, const vec& v);
    matrix& operator|=(matrix& A, const vec& v);

    // Vector operations
    vec add(const vec& v, const vec& w);
    vec operator+(const vec& v, const vec& w);
    void add_inplace(vec& v, const vec& w);
    vec& operator+=(vec& v, const vec& w);

    vec sub(const vec& v, const vec& w);
    vec operator-(const vec& v, const vec& w);
    void sub_inplace(vec& v, const vec& w);
    vec& operator-=(vec& v, const vec& w);

    vec scale(const vec& v, const double& sf);
    vec operator*(const vec& v, const double& sf);
    void scale_inplace(vec& v, const double& sf);
    vec& operator*=(vec& v, const double& sf);

    double dot(const vec& v, const vec& w);
    /// @brief Compute the dot product of 2 vectors
    double operator*(const vec& v, const vec& w);

    vec vecmult(const vec& v, const vec& w);
    /// @brief Compute the entry wise product of 2 vectors
    vec operator&(const vec& v, const vec& w);

    matrix dyadic(const vec& v, const vec& w);
    /// @brief Compute the dyadic product of 2 vectors
    matrix operator%(const vec& v, const vec& w);

    vec subvector(const vec& v, int start, int end); // End exclusive - as with iterators

    // --- Matrix operations
    // Basic
    matrix add(const matrix& A, const matrix& B);
    matrix operator+(const matrix& A, const matrix& B);
    void add_inplace(matrix& A, const matrix& B);
    matrix& operator+=(matrix& A, const matrix& B);

    matrix sub(const matrix& A, const matrix& B);
    matrix operator-(const matrix& A, const matrix& B);
    void sub_inplace(matrix& A, const matrix& B);
    matrix& operator-=(matrix& A, const matrix& B);

    matrix scale(const matrix& A, const double& sf);
    matrix operator*(const matrix& A, const double& sf);
    matrix& scale_inplace(matrix& A, const double& sf);
    matrix& operator*=(matrix& A, const double& sf);

    // Multiplications
    vec mult(const matrix& A, const vec& v);
    vec operator*(const matrix& A, const vec& v);

    matrix mult(const matrix& A, const matrix& B);
    matrix operator*(const matrix& A, const matrix& B);

    // Solving LSE's
    matrix inv(const matrix& A);
    matrix operator~(const matrix& A);

    double det(const matrix& A);

    vec solve(const matrix& A, const vec& v);
    vec operator/(const matrix& A, const vec& v);

    std::pair<double, matrix> gauss(const matrix& LHS, matrix&& RHS); // perform the gauss algorithm, getting det. and inv. together
    std::pair<double, matrix> gauss(const matrix& LHS, const vec& RHS);

    // TODO also include some special easier cases, like unit vectors
    /// @brief
    /// @return 
    matrix ShermanMorrison(const matrix& A, const vec& u, const vec& v);

    // Helper
    void swapRows(matrix& A, matrix& B, const int& i_A, const int& i_B);
} 

#endif // _LINALG