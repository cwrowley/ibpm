// n_gmres_vector
//
// Sample routine demonstrating the Newton-GMRES routine using a vector
// function.
//
// Kevin K. Chen
// Princeton University

#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <string>

#include <Eigen/Core>

#include "Gmres.h"
#include "Newton.h"
#include "NewtonArmijo.h"
#include "exceptions.h"
#include "n_gmres_vector.h"

int main(int argc, char* argv[]) {
    using std::cout;
    using std::cerr;
    using std::endl;
    using std::string;

    const string usage {
        "Usage: n_gmres_vector Newton | NewtonArmijo x0 x1"
    };

    if (argc != 4) {
        cerr << usage << endl;
        exit(1);
    }

    const string solver{argv[1]};
    Vec x{atof(argv[2]), atof(argv[3])};

    linear_solver_ns::Gmres<Vec> gmres(innerproduct, 2, 1.e-6);

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

    cout << "Initial guess: x = (" << x(0) << ", " << x(1) << ")" << endl;

    try {
        newton->solve(x);
    }
    catch (newton_ns::NewtonError& ne) {
        cerr << "Newton error: " << ne.what() << endl;
        exit(2);
    }
    catch (std::exception& e) {
        cerr << "Other exception: " << e.what() << endl;
        exit(3);
    }

    cout << "\nFinal solution: x = (" << x(0) << ", " << x(1) << ")." << endl;
}
