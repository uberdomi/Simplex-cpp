#ifndef _LINALG_BASIC
#define _LINALG_BASIC

#include "matrix.hpp"
#include <vector>

namespace la {
// --- Vector operations ---
template <typename NumType> void print(const std::vector<NumType> &v);
// Checks
template <typename NumType>
void check_dimensions(const std::vector<NumType> &v_left,
                      const std::vector<NumType> &v_right);
//   Operators
template <typename NumType>
std::vector<NumType> operator+(const std::vector<NumType> &v_left,
                               const std::vector<NumType> &v_right);
template <typename NumType>
std::vector<NumType> operator-(const std::vector<NumType> &v_left,
                               const std::vector<NumType> &v_right);
template <typename NumType>
std::vector<NumType> operator*(const std::vector<NumType> &v_left,
                               const std::vector<NumType> &v_right);
// Functions
template <typename NumType>
NumType dot(const std::vector<NumType> &v_left,
            const std::vector<NumType> &v_right);
template <typename NumType>
Matrix2D<NumType> dyadic(const std::vector<NumType> &v_left,
                         const std::vector<NumType> &v_right);
// --- Matrix operations ---
// Checks
template <typename NumType>
void check_dimensions(const Matrix2D<NumType> &m_left,
                      const Matrix2D<NumType> &m_right, uint dim);
// Dimensional operations
template <typename NumType> NumType sum(const Matrix2D<NumType> &A);
template <typename NumType>
inline std::vector<NumType> row_sum(const Matrix2D<NumType> &A);
template <typename NumType>
inline std::vector<NumType> col_sum(const Matrix2D<NumType> &A);
template <typename NumType>
Matrix2D<NumType> concat(const Matrix2D<NumType> &m_left,
                         const Matrix2D<NumType> &m_right, uint dim = 0);

// Matrix-vector
template <typename NumType>
std::vector<NumType> operator*(const std::vector<NumType> &v_left,
                               const std::vector<NumType> &v_right);

// Matrix-matrix
template <typename NumType>
Matrix2D<NumType> operator+(const Matrix2D<NumType> &m_left,
                            const Matrix2D<NumType> &m_right);

template <typename NumType>
Matrix2D<NumType> operator-(const Matrix2D<NumType> &m_left,
                            const Matrix2D<NumType> &m_right);

template <typename NumType>
Matrix2D<NumType> operator*(const Matrix2D<NumType> &m_left,
                            const Matrix2D<NumType> &m_right);

} // namespace la
#endif // _LINALG_BASIC