// continuationVector
//
// Sample routine demonstrating numerical continuation, using a vector
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
#include "continuationVector.h"

using namespace std::placeholders;

using continuation_ns::Continuation;
using eig_ns::AnalyzeEig;

int main(int argc, char* argv[]) {
    if (argc != 1) {
        std::cerr << "Too many input arguments." << std::endl;
        exit(1);
    }
    static_cast<void>(argv); // Turn off compiler warning about unused argv.

    constexpr int max_gmres_iter{2};
    constexpr double tol_gmres{1.e-6};
    constexpr double tol_newton{1.e-6};
    constexpr int max_newton_iter{25};
    constexpr double jacobian_dx{1.e-3};
    constexpr int jacobian_order{1};

    constexpr double lam_0{3.};
    constexpr double lam_f{-3.};
    constexpr double dlam_max{-1.};
    constexpr double dlam_min{-0.0625};
    const Vec x{-0.3, 0.2};

    linear_solver_ns::Gmres<Vec> gmres{innerProduct, max_gmres_iter, tol_gmres};
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
    AnalyzeEig<Vec> adjoint_analyze_eig{
        eigSolver,
        f_adjoint_lin,
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
        std::bind(&AnalyzeEig<Vec>::analyzeEig, &adjoint_analyze_eig, _1),
        param
    };

    continuation.solve();
}

void printVectorXcd(const VectorXcd& v, double lam) {
    using std::cout;
    using std::endl;

    cout << "Eigenvalues at lam = " << lam << ":" << endl;

    for (int j = {0}; j < v.size(); ++j) {
        cout << v(j).real();

        const double imag{v(j).imag()};
        if (imag < 0) {
            cout << " - " << -imag << "i" << endl;
        } else {
            cout << " + " << imag << "i" << endl;
        }
    }
}

void printEVectors(const vector<vector<Vec>>& eVectors, double lam) {
    using std::cout;
    using std::endl;

    assert(eVectors.size() == 2);
    unsigned n{static_cast<unsigned>(eVectors[0][0].size())};
    cout << "Eigenvectors at lam = " << lam << ":" << endl;

    for (unsigned i = {0}; i < eVectors.size(); ++i) {
        assert(eVectors[i].size() == 2);

        VectorXcd v{n};
        for (unsigned j = {0}; j < n; ++j) {
            v[j] = std::complex<double>{
                eVectors[i][0][j],
                eVectors[i][1][j]
            };
        }
        cout << v << "\n" << endl;
    }
}
