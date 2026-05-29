#ifndef _LINALG_BASIC
#define _LINALG_BASIC

#include "matrix.hpp"
#include <vector>

namespace la {
// --- Vector operations ---
void print(const std::vector<double> &v);
// Checks

void check_dimensions(const std::vector<double> &v_left,
                      const std::vector<double> &v_right);
//   Operators

std::vector<double> add(const std::vector<double> &v_left,
                        const std::vector<double> &v_right);

std::vector<double> sub(const std::vector<double> &v_left,
                        const std::vector<double> &v_right);

std::vector<double> mult(const std::vector<double> &v_left,
                         const std::vector<double> &v_right);
// Functions

double dot(const std::vector<double> &v_left,
           const std::vector<double> &v_right);

Matrix2D dyadic(const std::vector<double> &v_left,
                const std::vector<double> &v_right);
// --- Matrix operations ---
// Checks

void check_dimensions(const Matrix2D &m_left, const Matrix2D &m_right,
                      uint dim);
// Dimensional operations
double sum(const Matrix2D &A);

inline std::vector<double> row_sum(const Matrix2D &A);

inline std::vector<double> col_sum(const Matrix2D &A);

Matrix2D concat(const Matrix2D &m_left, const Matrix2D &m_right, uint dim = 0);

// Matrix-vector

std::vector<double> operator*(const Matrix2D &A, const std::vector<double> &v);

// Matrix-matrix

Matrix2D operator+(const Matrix2D &m_left, const Matrix2D &m_right);

Matrix2D operator-(const Matrix2D &m_left, const Matrix2D &m_right);

Matrix2D operator*(const Matrix2D &m_left, const Matrix2D &m_right);

} // namespace la
#endif // _LINALG_BASIC