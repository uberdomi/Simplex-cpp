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

    bool checkRectangular(const matrix& A);

    // Vector manipulations
    void scaleRow(vec& v, const double& sf);
    void subtractRow(vec& v, const vec& w, const double& sf=1.0);

    // Printing
    void print(const vec& v);
    void print(const matrix& A);

    // Conversions
    matrix vecAsMatrix(const vec& v);
    vec matrixAsVec(const matrix& A);
}

#endif // _UTIL