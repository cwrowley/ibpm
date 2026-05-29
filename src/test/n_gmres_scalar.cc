// n_gmres_scalar
//
// Sample routine demonstrating the Newton-Armijo-GMRES routine using a scalar
// function.
//
// Kevin K. Chen
// Princeton University

#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <string>

#include "Gmres.h"
#include "Newton.h"
#include "NewtonArmijo.h"
#include "exceptions.h"
#include "n_gmres_scalar.h"

int main(int argc, char* argv[]) {
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::string;

    const string usage{
        "Usage: n_gmres_scalar Newton | NewtonArmijo x"
    };

    if (argc != 3) {
        cerr << usage << endl;
        exit(1);
    }

    const string solver{argv[1]};
    Vec x{atof(argv[2])};

    linear_solver_ns::Gmres<Vec> gmres{ip, 1, 1.e-6};

    constexpr double tol{1.e-8};
    constexpr int max_iter{10};
    constexpr double jacobian_dx{1.e-3};

    std::unique_ptr<newton_ns::Newton<Vec>> newton;
    if (solver == "Newton") {
        newton.reset(
            new newton_ns::Newton<Vec>{f, gmres, norm, tol, max_iter,
                    jacobian_dx, true}
        );
    } else if (solver == "NewtonArmijo") {
        newton.reset(
            new newton_ns::NewtonArmijo<Vec>{f, gmres, norm, tol, max_iter,
                    jacobian_dx, 1, 10, 0.1, 0.5, 1.e-4, true}
        );
    } else {
        cerr << usage << endl;
        exit(2);
    }

    cout << "Initial guess: x = " << x << endl;

    try {
        newton->solve(x);
    }
    catch (newton_ns::NewtonError& ne) {
        cerr << "Newton error: " << ne.what() << endl;
        exit(3);
    }
    catch (std::exception& e) {
        cerr << "Other exception: " << e.what() << endl;
        exit(4);
    }

    cout << "\nFinal solution: x = " << x << endl;
}
