#ifndef _LINALG
#define _LINALG

#include <cstddef> // for std::size_t
#include <tuple>   // for std::tuple
#include <vector>

namespace la {

template <typename NumType> class Matrix2D {
private:
  std::size_t n_rows;
  std::size_t n_cols;

  std::size_t stride_rows;
  std::size_t stride_cols;

  std::vector<NumType>
      data; // All data stored contiguously, accesses defined by the strides

public:
  // --- Common constructors ---
  Matrix2D(std::size_t n_rows, std::size_t n_cols);

  Matrix2D(std::vector<std::vector<NumType>> inputs);

  // --- Stride operations ---
  inline NumType at(std::size_t row, std::size_t col);

  Matrix2D &t();

  // --- Getters ---
  std::size_t get_size();
  std::size_t get_n_rows();
  std::size_t get_n_cols();

  std::tuple<std::size_t, std::size_t> get_shape();
  std::tuple<std::size_t, std::size_t> get_strides();

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

#endif // _LINALG
