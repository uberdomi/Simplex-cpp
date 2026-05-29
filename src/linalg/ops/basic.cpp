#include "basic.hpp"

#include <algorithm>  // std::transform
#include <functional> // std::plus and other operations
#include <iostream>
#include <numeric> // std::accumulate
#include <stdexcept>

namespace la {
// --- Vector operations ---

template <typename NumType> void print(const std::vector<NumType> &v) {
  for (const auto &x : v) {
    std::cout << x << ", ";
  }
  std::cout << std::endl;
}

// Checks
template <typename NumType>
void check_dimensions(const std::vector<NumType> &v_left,
                      const std::vector<NumType> &v_right) {
  if (v_left.size() != v_right.size()) {
    throw std::runtime_error(
        "Incompatible vector dimensions: " + std::to_string(v_left.size()) +
        " vs " + std::to_string(v_right.size()));
  }
  if (v_left.empty()) {
    throw std::runtime_error("Vector is empty!");
  }
  return;
}
//   Operators

template <typename NumType>
std::vector<NumType> operator+(const std::vector<NumType> &v_left,
                               const std::vector<NumType> &v_right) {
  check_dimensions(v_left, v_right);
  std::vector<NumType> buffer = std::vector<NumType>(v_left.size(), NumType{});

  std::transform(v_left.begin(), v_left.end(), v_right.begin(), buffer.begin(),
                 [](NumType left, NumType right) { return left + right; });

  return buffer;
}

template <typename NumType>
std::vector<NumType> operator-(const std::vector<NumType> &v_left,
                               const std::vector<NumType> &v_right) {
  check_dimensions(v_left, v_right);

  std::vector<NumType> buffer = std::vector<NumType>(v_left.size(), NumType{});

  std::transform(v_left.begin(), v_left.end(), v_right.begin(), buffer.begin(),
                 [](NumType left, NumType right) { return left - right; });

  return buffer;
}
template <typename NumType>
std::vector<NumType> operator*(const std::vector<NumType> &v_left,
                               const std::vector<NumType> &v_right) {
  check_dimensions(v_left, v_right);

  std::vector<NumType> buffer = std::vector<NumType>(v_left.size(), NumType{});

  std::transform(v_left.begin(), v_left.end(), v_right.begin(), buffer.begin(),
                 [](NumType left, NumType right) { return left * right; });

  return buffer;
}

//  Functions

template <typename NumType>
NumType dot(const std::vector<NumType> &v_left,
            const std::vector<NumType> &v_right) {
  check_dimensions(v_left, v_right);

  size_t idx{0};

  NumType result = std::accumulate(v_left.begin(), v_left.end(), NumType{},
                                   [&v_right, &idx](NumType agg, NumType val) {
                                     return agg + val * v_right[idx++];
                                   });

  return result;
}

template <typename NumType>
Matrix2D<NumType> dyadic(const std::vector<NumType> &v_left,
                         const std::vector<NumType> &v_right) {
  if (v_left.empty()) {
    throw std::runtime_error("Left vector is empty!");
  }
  auto outer_vector = std::vector<std::vector<NumType>>{
      v_left.size(), std::vector<NumType>{v_right.size(), NumType{}}};
  for (size_t idx_left = 0; idx_left < v_left.size(); idx_left++) {
    for (size_t idx_right = 0; idx_right < v_right.size(); idx_right++) {
      outer_vector[idx_left][idx_right] = v_left[idx_left] * v_right[idx_right];
    }
  }

  return Matrix2D<NumType>(outer_vector);
}

//  --- Matrix operations ---
// Checks
template <typename NumType>
void check_dimensions(const Matrix2D<NumType> &m_left,
                      const Matrix2D<NumType> &m_right, uint dim) {
  switch (dim) {
  case 0: {
    if (m_left.get_n_rows() != m_right.get_n_rows()) {
      throw std::runtime_error("Incompatible matrix dimension 0: " +
                               std::to_string(m_left.get_n_rows()) + " vs " +
                               std::to_string(m_right.get_n_rows()));
    }
    break;
  }
  case 1: {
    if (m_left.get_n_cols() != m_right.get_n_cols()) {
      throw std::runtime_error("Incompatible matrix dimension 1: " +
                               std::to_string(m_left.get_n_cols()) + " vs " +
                               std::to_string(m_right.get_n_cols()));
    }
    break;
  }
  default: {
    if (m_left.get_n_cols() != m_right.get_n_cols() &&
        m_left.get_n_rows() != m_right.get_n_rows()) {
      throw std::runtime_error("Incompatible matrix dimensions: " +
                               std::to_string(m_left.get_shape()) + " vs " +
                               std::to_string(m_right.get_shape()));
    }
    break;
  }
  }
  return;
}

// Dimensional operations
template <typename NumType> NumType sum(const Matrix2D<NumType> &A, uint dim) {
  auto raw_data = A.get_raw_data;
  // Works cause whole data vector contiguous in memory - parts of it not
  // guaranteed! (vide A.transpose())
  NumType agg = std::accumulate(std::begin(raw_data), std::end(raw_data),
                                NumType{}, std::plus<NumType>{});
}

template <typename NumType>
inline std::vector<NumType> row_sum(const Matrix2D<NumType> &A) {
  std::vector<NumType> agg_vec =
      std::vector<NumType>{A.get_n_cols(), NumType{}};

  for (size_t idx_row = 0; idx_row < A.get_n_rows(); idx_row++) {
    for (size_t idx_col = 0; idx_col < A.get_n_cols(); idx_col++) {
      agg_vec[idx_row] += A(idx_row, idx_col);
    }
  }

  return agg_vec;
}

template <typename NumType>
inline std::vector<NumType> col_sum(const Matrix2D<NumType> &A) {
  std::vector<NumType> agg_vec =
      std::vector<NumType>{A.get_n_cols(), NumType{}};

  for (size_t idx_col = 0; idx_col < A.get_n_cols(); idx_col++) {
    for (size_t idx_row = 0; idx_row < A.get_n_rows(); idx_row++) {
      agg_vec[idx_col] += A(idx_row, idx_col);
    }
  }

  return agg_vec;
}

} // namespace la