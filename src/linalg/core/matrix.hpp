#ifndef _LINALG_MATRIX
#define _LINALG_MATRIX

#include <cstddef> // for std::size_t
#include <memory>  // for pointers
#include <utility> // for std::pair
#include <vector>

namespace la {

template <typename NumType> class Matrix2D {
private:
  // Everything immutable, but new objects can be created out of it
  std::size_t n_rows_;
  std::size_t n_cols_;

  std::size_t stride_rows_;
  std::size_t stride_cols_;

  std::size_t size_;

  std::shared_ptr<const NumType[]>
      data_; // All data stored contiguously, accesses defined by the strides

public:
  // --- Public constructors ---
  Matrix2D(std::size_t n_rows, std::size_t n_cols);

  Matrix2D(std::vector<std::vector<NumType>> inputs);

private:
  // Constructing a matrix from the same underlying data *without copying*
  Matrix2D(std::shared_ptr<const NumType[]> shared_data, std::size_t n_rows,
           std::size_t n_cols);

public:
  // --- Stride operations ---

  // Raw access - no bounds checking
  inline NumType operator()(std::size_t row, std::size_t col);

  // Civilized access with bound checking
  inline NumType at(std::size_t row, std::size_t col);

  Matrix2D transpose();

  // --- Getters ---
  std::size_t get_size();
  std::size_t get_n_rows();
  std::size_t get_n_cols();

  std::pair<std::size_t, std::size_t> get_shape();
  std::pair<std::size_t, std::size_t> get_strides();

  std::shared_ptr<const NumType[]> get_raw_data();

  // --- Further functionalities ---
  void print();
};

// --- Common matrices ---

Matrix2D<double> zeros(std::size_t n_rows, std::size_t n_cols);
Matrix2D<double> zeros(std::size_t length);

Matrix2D<double> ones(std::size_t n_rows, std::size_t n_cols);
Matrix2D<double> ones(std::size_t length);

Matrix2D<double> eye(std::size_t length);

} // namespace la

#endif // _LINALG_MATRIX