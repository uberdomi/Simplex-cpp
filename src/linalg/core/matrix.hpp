#ifndef _LINALG
#define _LINALG

#include <vector>

namespace la {

template <typename NumType> class Matrix2D {
private:
  size_t n_rows;
  size_t n_cols;

  size_t stride_rows;
  size_t stride_cols;

  std::vector<NumType>
      data; // All data stored contiguously, accesses defined by the strides

public:
  // --- Common constructors ---
  Matrix2D(size_t n_rows, size_t n_cols);

  Matrix2D(std::vector<std::vector<NumType>> inputs);

  // --- Stride operations ---
  NumType at(std::size_t row, std::size_t col);

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

Matrix2D<double> zeros(size_t n_rows, size_t n_cols);
Matrix2D<double> zeros(size_t length);

Matrix2D<double> ones(size_t n_rows, size_t n_cols);
Matrix2D<double> ones(size_t length);

Matrix2D<double> eye(size_t length);

} // namespace la

#endif // _LINALG
