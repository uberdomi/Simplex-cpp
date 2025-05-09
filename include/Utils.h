#ifndef _UTIL
#define _UTIL

#include <vector>

namespace util {
    // Helper structures
    using vec = std::vector<double>;
    using matrix = std::vector<vec>;

    // Comparisons
    bool isEqual(const double& x, const double& y, const double& tol=1e-6);
    bool isEqual(const vec& v, const vec& w, const double& tol=1e-6);
    bool isEqual(const matrix& A, const matrix& B, const double& tol=1e-6);

    bool isRectangular(const matrix& A);

    // Vector manipulations
    void scaleRow(vec& v, const double& sf);
    void subtractRow(vec& v, const vec& w, const double& sf=1.0);

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

    void append(matrix& A, matrix&& B);
    void append(matrix& A, vec&& v);
    void append(vec& v, vec&& w);

    // Vector operations
    vec add(const vec& v, const vec& w);
    vec sub(const vec& v, const vec& w);
    vec scale(const vec& v, const double& sf);

    double dot(const vec& v, const vec& w);
    matrix dyadic(const vec& v, const vec& w);

    vec subvector(const vec& v, int start, int end); // End exclusive - as with iterators
    
}

#endif // _UTIL