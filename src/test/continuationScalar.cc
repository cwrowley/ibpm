// continuationScalar
//
// Sample routine demonstrating numerical continuation, using a scalar
// function.
//
// Kevin K. Chen
// Princeton University

#include <cstdlib>
#include <functional>
#include <iostream>

#include <Eigen>
#include <Core>

#include "AnalyzeEig.h"
#include "Continuation.h"
#include "Newton.h"
#include "Gmres.h"
#include "continuationScalar.h"

using std::cerr;
using std::endl;
using namespace std::placeholders;

using continuation_ns::Continuation;
using eig_ns::AnalyzeEig;

int main(int argc, char* argv[]) {
    if (argc != 6) {
        cerr << "Usage: continuationScalar IC lam_0 lam_f dlam_max dlam_min"
             << endl;
        exit(1);
    }

    constexpr int max_gmres_iter{1};
    constexpr double tol_gmres{1.e-6};
    constexpr double tol_newton{1.e-6};
    constexpr int max_newton_iter{10};
    constexpr double jacobian_dx{1.e-3};
    constexpr int jacobian_order{1};

    const Vec x{atof(argv[1])};
    const double lam_0{atof(argv[2])};
    const double lam_f{atof(argv[3])};
    const double dlam_max{atof(argv[4])};
    const double dlam_min{atof(argv[5])};

    linear_solver_ns::Gmres<Vec> gmres{ip, max_gmres_iter, tol_gmres};
    newton_ns::Newton<Vec> newton{
        gmres,
        norm,
        tol_newton,
        max_newton_iter,
        jacobian_dx,
        jacobian_order,
        true
    };

    Continuation<Vec>::Param param{"lambda", lam_0, lam_f, dlam_max, dlam_min};
    AnalyzeEig<Vec> analyze_eig{
        eigSolver,
        f_lin,
        eig_ns::OperatorType::continuous,
        printVectorXcd,
        printEVectors,
        x,
        param.lam
    };
    Continuation<Vec> continuation{
        f,
        newton,
        x,
        print,
        std::bind(&AnalyzeEig<Vec>::analyzeEig, &analyze_eig, _1),
        std::bind(&AnalyzeEig<Vec>::analyzeEig, &analyze_eig, _1),
        param
    };

    continuation.solve();
}

void printVectorXcd(const VectorXcd& v, double lam) {
    cout << "Eigenvalues at lam = " << lam << ":" << endl;

    for (int j = {0}; j < v.size(); ++j) {
        cout << v(j).real();

        double imag = v(j).imag();
        if (imag < 0) {
            cout << " - " << -imag << "i" << endl;
        } else {
            cout << " + " << imag << "i" << endl;
        }
    }
}

void printEVectors(const vector<vector<Vec>>& evectors, double lam) {
    assert(evectors.size() == 1);
    assert(evectors[0].size() == 2);

    cout << "Eigenvector at lam = " << lam << ":" << endl;
    cout << evectors[0][0] << " + " << evectors[0][1] << "i" << endl;
}
