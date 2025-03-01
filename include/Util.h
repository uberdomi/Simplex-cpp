#include <vector>

namespace util {
    using vec = std::vector<double>;
    using matrix = std::vector<vec>;

    bool isEqual(const double&, const double&, const double& tol=1e-6);
    bool isEqual(const vec&, const vec&, const double& tol=1e-6);
    bool isEqual(const matrix&, const matrix&, const double& tol=1e-6);

    void scaleRow(vec&, double);
    void subtractRow(vec&, vec&, double sf=1.0);

    void print(const vec&);
    void print(const matrix&);
}