#ifndef N_GMRES_VECTOR_H
#define N_GMRES_VECTOR_H

#include <Eigen/Core>

using Vec = Eigen::Vector2d;

inline double innerproduct(const Vec& x, const Vec& y) {
    return y.dot(x);
}

inline double norm(const Vec& x) {
    return x.norm();
}

// Sample nonlinear function.
inline Vec f(const Vec& x) {
    static int num_evals = 0;

    Vec y{x(0) * x(0) - 9 + x(1), 3 * x(0) + x(1) - 5};

    std::cout << "\t[" << ++num_evals << "] f(" << x(0) << ", " << x(1)
              << ") = " << "(" << y(0) << ", " << y(1) << ")" << std::endl;

    return y;
}

#endif // N_GMRES_VECTOR_H
