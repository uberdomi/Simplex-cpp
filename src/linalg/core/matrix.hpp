#ifndef _LINALG
#define _LINALG

#include <vector>

namespace la {

template <typename Type>
class Matrix2D
{
private:
    size_t n_0;
    size_t n_1;

    size_t stride_0;
    size_t stride_1;

public:
    std::vector<Type> data; // All data stored contiguously, accesses defined by the strides
    std::tuple<size_t, size_t> shape;
    std::tuple<size_t, size_t> strides;
    size_t size = n_0 * n_1;

    Matrix2D(size_t n_rows, size_t n_cols) : n_0(n_rows), n_1(n_cols), stride_0(n_1), stride_1(1), data(n_rows * n_cols, Type()), shape(n_0, n_1), strides(stride_0, stride_1)
    {}

    
};

} // la





#endif // _LINALG
