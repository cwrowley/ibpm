#ifndef N_GMRES_SCALAR_H
#define N_GMRES_SCALAR_H

#include <cmath>
#include <iostream>

using Vec = double;

// Sample nonlinear function.
inline Vec f(const Vec& x) {
    static int numEvals = 0;
    Vec y{std::atan(x)};
    std::cout << "\t[" << ++numEvals << "] f(" << x << ") = " << y << std::endl;

    return y;
}

inline double ip(const Vec& x, const Vec& y){
    return x * y;
}

inline double norm(const Vec& x) {
    return sqrt(ip(x, x));
}

#endif // N_GMRES_SCALAR_H
