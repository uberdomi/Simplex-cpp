#include "matrix.hpp"

#include <iostream>
#include <stdexcept>

// --- Constructors ---

namespace la {

template <typename NumType>
Matrix2D<NumType>::Matrix2D(std::size_t n_rows, std::size_t n_cols)
    : n_rows(n_rows), n_cols(n_cols), stride_rows(n_cols), stride_cols(1),
      data(n_rows * n_cols, NumType()) {}

template <typename NumType>
Matrix2D<NumType>::Matrix2D(std::vector<std::vector<NumType>> inputs) {
  if (inputs.empty()) {
    throw std::runtime_error("Empty input matrix data!");
  }
  n_rows = inputs.size();

  if (inputs[0].empty()) {
    throw std::runtime_error("Empty input matrix row!");
  }

  n_cols = inputs[0].size();

  stride_rows = n_cols;
  stride_cols = 1;

  data = std::vector(n_rows * n_cols, NumType());

  std::size_t offset = 0;

  for (const auto &row : inputs) {
    if (row.size() != n_cols) {
      throw std::runtime_error("Row size mismatch : current " +
                               std::to_string(row.size()) + " vs wanted " +
                               std::to_string(n_cols));
    }

    std::copy(row.begin(), row.end(), data.begin() + offset);

    offset += n_cols;
  }
}

// --- Stride operations ---

template <typename NumType>
inline NumType Matrix2D<NumType>::at(std::size_t row, std::size_t col) {
  //   Automatic bounds checking by the in-build function
  return data.at(row * stride_rows + col * stride_cols);
}

template <typename NumType> Matrix2D<NumType> &Matrix2D<NumType>::t() {
  //   Automatic bounds checking by the in-build function
  std::swap(stride_rows, stride_cols);

  return *this;
}

// --- Getters ---

template <typename NumType> std::size_t Matrix2D<NumType>::get_size() {
  return n_rows * n_cols;
}

template <typename NumType> std::size_t Matrix2D<NumType>::get_n_rows() {
  return n_rows;
}

template <typename NumType> std::size_t Matrix2D<NumType>::get_n_cols() {
  return n_cols;
}

template <typename NumType>
std::tuple<std::size_t, std::size_t> Matrix2D<NumType>::get_shape() {
  return {n_rows, n_cols};
}

template <typename NumType>
std::tuple<std::size_t, std::size_t> Matrix2D<NumType>::get_strides() {
  return {stride_rows, stride_cols};
}

// --- Further functionalities ---

template <typename NumType> void Matrix2D<NumType>::print() {
  std::cout << "--- matrix2d (" << n_rows << "," << n_cols << ")" << std::endl;
  for (std::size_t row = 0; row < n_rows; row++) {
    std::cout << "[ ";
    for (std::size_t col = 0; col < n_cols; col++) {
      std::cout << this->at(row, col) << ", ";
    }
    std::cout << "]" << std::endl;
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