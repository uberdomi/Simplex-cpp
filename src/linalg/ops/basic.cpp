#include "basic.hpp"

#include <algorithm> // std::transform
#include <iostream>
#include <numeric> // std::accumulate
#include <span>    // for std::span
#include <stdexcept>

namespace la {
// --- Vector operations ---

void print(const std::vector<double> &v) {
  for (const auto &x : v) {
    std::cout << x << ", ";
  }
  std::cout << std::endl;
}

// Checks

void check_dimensions(const std::vector<double> &v_left,
                      const std::vector<double> &v_right) {
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

std::vector<double> add(const std::vector<double> &v_left,
                        const std::vector<double> &v_right) {
  check_dimensions(v_left, v_right);
  std::vector<double> buffer = std::vector<double>(v_left.size(), 0.0);

  std::transform(v_left.begin(), v_left.end(), v_right.begin(), buffer.begin(),
                 [](double left, double right) { return left + right; });

  return buffer;
}

std::vector<double> sub(const std::vector<double> &v_left,
                        const std::vector<double> &v_right) {
  check_dimensions(v_left, v_right);

  std::vector<double> buffer = std::vector<double>(v_left.size(), 0.0);

  std::transform(v_left.begin(), v_left.end(), v_right.begin(), buffer.begin(),
                 [](double left, double right) { return left - right; });

  return buffer;
}

std::vector<double> mult(const std::vector<double> &v_left,
                         const std::vector<double> &v_right) {
  check_dimensions(v_left, v_right);

  std::vector<double> buffer = std::vector<double>(v_left.size(), 0.0);

  std::transform(v_left.begin(), v_left.end(), v_right.begin(), buffer.begin(),
                 [](double left, double right) { return left * right; });

  return buffer;
}

//  Functions

double dot(const std::vector<double> &v_left,
           const std::vector<double> &v_right) {
  check_dimensions(v_left, v_right);

  size_t idx{0};

  double result = std::accumulate(v_left.begin(), v_left.end(), 0.0,
                                  [&v_right, &idx](double agg, double val) {
                                    return agg + val * v_right[idx++];
                                  });

  return result;
}

Matrix2D dyadic(const std::vector<double> &v_left,
                const std::vector<double> &v_right) {
  if (v_left.empty()) {
    throw std::runtime_error("Left vector is empty!");
  }
  std::vector<std::vector<double>> outer_vector =
      std::vector<std::vector<double>>{
          v_left.size(), std::vector<double>(v_right.size(), 0.0)};
  for (size_t idx_left = 0; idx_left < v_left.size(); idx_left++) {
    for (size_t idx_right = 0; idx_right < v_right.size(); idx_right++) {
      outer_vector[idx_left][idx_right] = v_left[idx_left] * v_right[idx_right];
    }
  }

  return Matrix2D(outer_vector);
}

//  --- Matrix operations ---
// Checks

void check_dimensions(const Matrix2D &m_left, const Matrix2D &m_right,
                      uint dim) {
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
      throw std::runtime_error(
          "Incompatible matrix dimensions: " + m_left.shape_str() + " vs " +
          m_right.shape_str());
    }
    break;
  }
  }
  return;
}

// Dimensional operations
double sum(const Matrix2D &A, uint dim) {
  std::shared_ptr<const double[]> raw_data = A.get_raw_data();
  size_t data_size = A.get_size();

  std::span<const double> data_view(raw_data.get(), data_size);

  // Works cause whole data vector contiguous in memory - parts of it not
  // guaranteed! (vide A.transpose())
  double agg = std::accumulate(std::begin(data_view), std::end(data_view), 0.0);

  return agg;
}

inline std::vector<double> row_sum(const Matrix2D &A) {
  std::vector<double> agg_vec = std::vector<double>(A.get_n_cols(), 0.0);

  for (size_t idx_row = 0; idx_row < A.get_n_rows(); idx_row++) {
    for (size_t idx_col = 0; idx_col < A.get_n_cols(); idx_col++) {
      agg_vec[idx_row] += A(idx_row, idx_col);
    }
  }

  return agg_vec;
}

inline std::vector<double> col_sum(const Matrix2D &A) {
  std::vector<double> agg_vec = std::vector<double>(A.get_n_cols(), 0.0);

  for (size_t idx_col = 0; idx_col < A.get_n_cols(); idx_col++) {
    for (size_t idx_row = 0; idx_row < A.get_n_rows(); idx_row++) {
      agg_vec[idx_col] += A(idx_row, idx_col);
    }
  }

  return agg_vec;
}

} // namespace la