#include "matrix.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

// --- Constructors ---

namespace la {

// Default constructor
template <typename NumType>
Matrix2D<NumType>::Matrix2D(std::size_t n_rows, std::size_t n_cols)
    : n_rows_(n_rows), n_cols_(n_cols), stride_rows_(n_cols), stride_cols_(1),
      size_(n_rows * n_cols),
      data_(std::make_shared<NumType[]>(n_rows * n_cols)) {}

// Matrix from a 2D vector
template <typename NumType>
Matrix2D<NumType>::Matrix2D(std::vector<std::vector<NumType>> inputs) {
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
  std::shared_ptr<NumType[]> buffer = std::make_shared<NumType[]>(size_);

  for (size_t i; i < n_rows_; i++) {
    for (size_t j; j < n_cols_; j++) {
      buffer[i * stride_rows_ + j * stride_cols_] = inputs[i][j];
    }
  }

  data_ = std::move(buffer);
}

// Private constructor from a data array
template <typename NumType>
Matrix2D<NumType>::Matrix2D(std::shared_ptr<const NumType[]> shared_data,
                            std::size_t n_rows, std::size_t n_cols)
    : n_rows_(n_rows), n_cols_(n_cols), stride_rows_(n_cols), stride_cols_(1),
      size_(n_rows * n_cols), data_(std::move(shared_data)) {}

// --- Stride operations ---

template <typename NumType>
inline NumType Matrix2D<NumType>::operator()(std::size_t row, std::size_t col) {
  return data_[row * stride_rows_ + col * stride_cols_];
}

template <typename NumType>
inline NumType Matrix2D<NumType>::at(std::size_t row, std::size_t col) {
  //   With bounds checking
  size_t index = row * stride_rows_ + col * stride_cols_;
  if (index >= size_ || index < 0) {
    throw std::out_of_range("Invalid access index: " + std::to_string(index) +
                            " not in (0, " + std::to_string(size_) + ")");
  }
  return data_[index];
}

template <typename NumType> Matrix2D<NumType> Matrix2D<NumType>::transpose() {
  // Swap the rows and cols dimensions, share the underlying data
  return Matrix2D<NumType>(data_, n_cols_, n_rows_);
}

// --- Getters ---

template <typename NumType> std::size_t Matrix2D<NumType>::get_size() {
  return size_;
}

template <typename NumType> std::size_t Matrix2D<NumType>::get_n_rows() {
  return n_rows_;
}

template <typename NumType> std::size_t Matrix2D<NumType>::get_n_cols() {
  return n_cols_;
}

template <typename NumType>
std::pair<std::size_t, std::size_t> Matrix2D<NumType>::get_shape() {
  return {n_rows_, n_cols_};
}

template <typename NumType>
std::pair<std::size_t, std::size_t> Matrix2D<NumType>::get_strides() {
  return {stride_rows_, stride_cols_};
}

template <typename NumType>
std::shared_ptr<const NumType[]> Matrix2D<NumType>::get_raw_data() {
  return data_;
}

// --- Further functionalities ---

template <typename NumType> void Matrix2D<NumType>::print() {
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

// Needed to be used for the correct return types
template class Matrix2D<double>;

// --- Common matrices ---

// Zeros
Matrix2D<double> zeros(std::size_t n_rows, std::size_t n_cols) {
  std::vector<std::vector<double>> inputs(n_rows,
                                          std::vector<double>(n_cols, 0.0));

  return Matrix2D(inputs);
}

Matrix2D<double> zeros(std::size_t length) { return zeros(length, length); }

// Ones
Matrix2D<double> ones(std::size_t n_rows, std::size_t n_cols) {
  std::vector<std::vector<double>> inputs(n_rows,
                                          std::vector<double>(n_cols, 1.0));

  return Matrix2D(inputs);
}
Matrix2D<double> ones(std::size_t length) { return ones(length, length); }

// Identity
Matrix2D<double> eye(std::size_t length) {
  std::vector<std::vector<double>> inputs(length,
                                          std::vector<double>(length, 0.0));

  for (std::size_t idx = 0; idx < length; idx++) {
    inputs[idx][idx] = 1.0;
  }

  return Matrix2D(inputs);
}

} // namespace la