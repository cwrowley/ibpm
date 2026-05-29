#ifndef CONTINUATION_SCALAR_H
#define CONTINUATION_SCALAR_H

#include <cassert>
#include <iostream>
#include <vector>

#include <Eigen>
#include <Core>

#include "ArnoldiEig.h"
#include "templateAliases.h"

using std::cout;
using std::endl;
using std::vector;

using Eigen::VectorXcd;

using Vec = double;

inline Vec f(const Vec& x, double lam) {
    static int numEvals{0};

    Vec fx{-x * x * x + x * lam}; // Supercritical pitchfork.

    cout << "\t[" << ++numEvals << "] f(" << x << ", " << lam << ") = " << fx
         << endl;

    return fx;
}

// Linearized version of f about an equilibrium x_eq.
inline Vec f_lin(const Vec& x, const Vec& x_eq, double lam) {
    return (-3 * x_eq * x_eq + lam) * x;
}

inline double ip(const Vec& x, const Vec& y) {
    return x * y;
}

inline double norm(const Vec& x) {
    return sqrt(ip(x, x));
}

inline void eigSolver(
    alias::Function<Vec> f,
    VectorXcd& eValues,
    vector<vector<Vec>>& evectors
) {
    eig_ns::ArnoldiEig<Vec> ae{f, ip, 1., 1, eig_ns::SortByType::real, true};
    ae.eig(eValues, evectors);
}

inline void print(const Vec& x, double lam) {
    cout << "\nSolution found: x = " << x << " at lambda = " << lam << "."
         << endl;
}

void printVectorXcd(const VectorXcd& v, double lam);
void printEVectors(const vector<vector<Vec>>& evectors, double lam);

#endif // CONTINUATION_SCALAR_H
