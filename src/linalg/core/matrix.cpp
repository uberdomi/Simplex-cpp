#include "matrix.hpp"

#include <iostream>
#include <stdexcept>

// --- Constructors ---

namespace la {

template <typename NumType>
Matrix2D<NumType>::Matrix2D(size_t n_rows, size_t n_cols)
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

  size_t offset = 0;

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
NumType Matrix2D<NumType>::at(std::size_t row, std::size_t col) {
  //   Automatic bounds checking by the in-build function
  return data.at(row * stride_rows + col * stride_cols);
}

template <typename NumType> Matrix2D<NumType> &Matrix2D<NumType>::t() {
  //   Automatic bounds checking by the in-build function
  std::swap(stride_rows, stride_cols);

  return this;
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
  std::cout << "{" << std::endl;
  for (const auto &row : data) {
    std::cout << "[ ";
    for (NumType val : row) {
      std::cout << val << ", ";
    }
    std::cout << "]" << std::endl;
  }
  std::cout << "}" << std::endl;
}

// --- Common matrices ---

// Zeros
Matrix2D<double> zeros(size_t n_rows, size_t n_cols) {
  std::vector<std::vector<double>> inputs(n_rows,
                                          std::vector<double>(n_cols, 0.0));

  return Matrix2D(inputs);
}

Matrix2D<double> zeros(size_t length) { return zeros(length, length); }

// Ones
Matrix2D<double> ones(size_t n_rows, size_t n_cols) {
  std::vector<std::vector<double>> inputs(n_rows,
                                          std::vector<double>(n_cols, 1.0));

  return Matrix2D(inputs);
}
Matrix2D<double> ones(size_t length) { return ones(length, length); }

// Identity
Matrix2D<double> eye(size_t length) {
  std::vector<std::vector<double>> inputs(length,
                                          std::vector<double>(length, 0.0));

  for (size_t idx = 0; idx < length; idx++) {
    inputs[idx][idx] = 1.0;
  }

  return Matrix2D(inputs);
}

} // namespace la