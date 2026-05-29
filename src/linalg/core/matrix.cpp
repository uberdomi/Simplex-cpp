#include "matrix.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

// --- Constructors ---

namespace la {

// Default constructor
Matrix2D::Matrix2D(std::size_t n_rows, std::size_t n_cols)
    : n_rows_(n_rows), n_cols_(n_cols), stride_rows_(n_cols), stride_cols_(1),
      size_(n_rows * n_cols),
      data_(std::make_shared<double[]>(n_rows * n_cols)) {}

// Matrix from a 2D vector
Matrix2D::Matrix2D(std::vector<std::vector<double>> inputs) {
  // Size checks
  if (inputs.empty()) {
    throw std::runtime_error("Empty input matrix data!");
  }

  if (inputs[0].empty()) {
    throw std::runtime_error("Empty input matrix row!");
  }

  n_rows_ = inputs.size();

  n_cols_ = inputs[0].size();

  stride_rows_ = n_cols_;
  stride_cols_ = 1;
  size_ = n_rows_ * n_cols_;

  // Check loop
  for (const auto &row_vec : inputs) {
    if (row_vec.size() != n_cols_) {
      throw std::runtime_error("Row size mismatch : current " +
                               std::to_string(row_vec.size()) + " vs wanted " +
                               std::to_string(n_cols_));
    }
  }

  // Populate the underlying data
  std::shared_ptr<double[]> buffer_ptr = std::make_shared<double[]>(size_);

  // auto buffer = buffer_ptr.get();

  for (std::size_t row_idx = 0; row_idx < n_rows_; row_idx++) {
    for (std::size_t col_idx = 0; col_idx < n_cols_; col_idx++) {
      buffer_ptr[row_idx * stride_rows_ + col_idx * stride_cols_] =
          inputs[row_idx][col_idx];
    }
  }

  data_ = buffer_ptr;
}

// Private constructor from a data array
Matrix2D::Matrix2D(std::shared_ptr<const double[]> shared_data,
                   std::size_t n_rows, std::size_t n_cols)
    : n_rows_(n_rows), n_cols_(n_cols), stride_rows_(n_cols), stride_cols_(1),
      size_(n_rows * n_cols), data_(std::move(shared_data)) {}

// --- Stride operations ---
inline double Matrix2D::operator()(std::size_t row, std::size_t col) const {
  return data_[row * stride_rows_ + col * stride_cols_];
}

inline double Matrix2D::at(std::size_t row, std::size_t col) const {
  //   With bounds checking
  size_t index = row * stride_rows_ + col * stride_cols_;
  if (index >= size_ || index < 0) {
    throw std::out_of_range("Invalid access index: " + std::to_string(index) +
                            " not in (0, " + std::to_string(size_) + ")");
  }
  return data_[index];
}
Matrix2D Matrix2D::transpose() {
  // Swap the rows and cols dimensions, share the underlying data
  return Matrix2D(data_, n_cols_, n_rows_);
}

// --- Getters ---
std::size_t Matrix2D::get_size() const { return size_; }
std::size_t Matrix2D::get_n_rows() const { return n_rows_; }
std::size_t Matrix2D::get_n_cols() const { return n_cols_; }

std::pair<std::size_t, std::size_t> Matrix2D::get_shape() const {
  return {n_rows_, n_cols_};
}

std::pair<std::size_t, std::size_t> Matrix2D::get_strides() const {
  return {stride_rows_, stride_cols_};
}

std::shared_ptr<const double[]> Matrix2D::get_raw_data() const { return data_; }

// --- Further functionalities ---
void Matrix2D::print() const {
  std::cout << "--- matrix2d (" << n_rows_ << "," << n_cols_ << ")"
            << std::endl;
  for (std::size_t row = 0; row < n_rows_; row++) {
    std::cout << "[ ";
    for (std::size_t col = 0; col < n_cols_; col++) {
      std::cout << this->at(row, col) << ", ";
    }
    std::cout << "]\n";
  }
  std::cout << "-------------" << std::endl;
}

std::string Matrix2D::shape_str() const {
  return "(" + std::to_string(n_rows_) + "," + std::to_string(n_cols_) + ")";
}

// --- Common matrices ---

// Zeros
Matrix2D zeros(std::size_t n_rows, std::size_t n_cols) {
  std::vector<std::vector<double>> inputs(n_rows,
                                          std::vector<double>(n_cols, 0.0));

  return Matrix2D(inputs);
}

Matrix2D zeros(std::size_t length) { return zeros(length, length); }

// Ones
Matrix2D ones(std::size_t n_rows, std::size_t n_cols) {
  std::vector<std::vector<double>> inputs(n_rows,
                                          std::vector<double>(n_cols, 1.0));

  return Matrix2D(inputs);
}
Matrix2D ones(std::size_t length) { return ones(length, length); }

// Identity
Matrix2D eye(std::size_t length) {
  std::vector<std::vector<double>> inputs(length,
                                          std::vector<double>(length, 0.0));

  for (std::size_t idx = 0; idx < length; idx++) {
    inputs[idx][idx] = 1.0;
  }

  return Matrix2D(inputs);
}

} // namespace la