#ifndef CONTINUATION_VECTOR_H
#define CONTINUATION_VECTOR_H

#include <functional>
#include <iostream>
#include <vector>

#include <Eigen>
#include <Core>
#include "ArnoldiEig.h"
#include "random.h"
#include "templateAliases.h"

using std::function;
using std::vector;

using Eigen::VectorXcd;

using Vec = Eigen::Vector2d;

inline Vec f(const Vec& x, double lam) {
    static int numEvals{0};

    // Supercritical Hopf bifurcation.
    double z{lam - x(0) * x(0) - x(1) * x(1)};
    Vec fx{-x(1) + x(0) * z, x(0) + x(1) * z};

    std::cout << "\t[" << ++numEvals << "] f((" << x.transpose() << "), " << lam
              << ") = (" << fx.transpose() << ")" << std::endl;

    return fx;
}

// Linearized version of f about an equilibrium x_eq.
inline Vec f_lin(const Vec& x, const Vec& x_eq, double lam) {
    Eigen::Matrix2d A{};
    A << lam, -1, 1., lam;
    static_cast<void>(x_eq); // Unused.
    return A * x;
}

// Adjoint linearized version of f about an equilibrium x_eq.
inline Vec f_adjoint_lin(const Vec& x, const Vec& x_eq, double lam) {
    Eigen::Matrix2d A{};
    A << lam, 1, -1., lam;
    static_cast<void>(x_eq); // Unused.
    return A * x;
}

inline double norm(const Vec& x) {
    return x.norm();
}

inline double innerProduct(const Vec& x, const Vec& y) {
    return y.dot(x);
}

inline void eigSolver(
    alias::Function<Vec> f,
    VectorXcd& eValues,
    vector<vector<Vec>>& eVectors
) {
    constexpr int n{2};
    Vec initial(n);
    for (int i = 0; i < n; ++i) {
        initial[i] = randDouble(-0.5, 0.5);
    }

    eig_ns::ArnoldiEig<Vec> ae{
        f,
        innerProduct,
        initial,
        n,
        eig_ns::SortByType::real,
        true
    };
    ae.eig(eValues, eVectors);
}

inline void print(const Vec& x, double lam) {
    std::cout << "\nSolution found: x = (" << x.transpose() << ") at lambda = "
              << lam << "." << std::endl;
}

void printVectorXcd(const VectorXcd& v, double lam);
void printEVectors(const vector<vector<Vec>>& eVectors, double lam);

#endif // CONTINUATION_VECTOR_H
