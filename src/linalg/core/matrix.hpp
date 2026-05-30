#ifndef _LINALG_MATRIX
#define _LINALG_MATRIX

#include <cstddef> // for std::size_t
#include <memory>  // for pointers
#include <string>
#include <utility> // for std::pair
#include <vector>

namespace la {

class Matrix2D {
private:
  // Everything immutable, but new objects can be created out of it
  std::size_t n_rows_;
  std::size_t n_cols_;

  std::size_t stride_rows_;
  std::size_t stride_cols_;

  std::size_t size_;

  std::shared_ptr<const double[]>
      data_; // All data stored contiguously, accesses defined by the strides

public:
  // --- Public constructors ---
  Matrix2D(std::size_t n_rows, std::size_t n_cols);

  Matrix2D(std::vector<std::vector<double>> inputs);

  // Constructing a matrix from the same underlying data *without copying*
  Matrix2D(std::shared_ptr<const double[]> shared_data, std::size_t n_rows,
           std::size_t n_cols);

  // --- Stride operations ---

  // Raw access - no bounds checking
  double operator()(std::size_t row, std::size_t col) const;

  // Civilized access with bound checking
  double at(std::size_t row, std::size_t col) const;

  Matrix2D transpose();

  // --- Getters ---
  std::size_t get_size() const;
  std::size_t get_n_rows() const;
  std::size_t get_n_cols() const;

  std::pair<std::size_t, std::size_t> get_shape() const;
  std::pair<std::size_t, std::size_t> get_strides() const;

  std::shared_ptr<const double[]> get_raw_data() const;

  // --- Further functionalities ---
  void print() const;
  std::string shape_str() const;
};

// --- Common matrices ---

Matrix2D zeros(std::size_t n_rows, std::size_t n_cols);
Matrix2D zeros(std::size_t length);

Matrix2D ones(std::size_t n_rows, std::size_t n_cols);
Matrix2D ones(std::size_t length);

Matrix2D eye(std::size_t length);

} // namespace la

#endif // _LINALG_MATRIX