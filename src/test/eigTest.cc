// eigTest - Test eigenvalue and eigenvector computation using Arnoldi.

#include "eigTest.h"

#include <complex>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen>
#include <Core>
#include <Eigenvalues>

#include "ArnoldiEig.h"
#include "random.h"

int main(int argc, char* argv[]) {
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::string;
    using std::vector;
    using namespace std::placeholders;

    using Eigen::MatrixXcd;
    using Eigen::VectorXcd;
    using Eigen::VectorXd;

    using eig_ns::SortByType;

    const string usage{"Usage: eigTest <size> nEVectors real | abs"};

    if (argc != 4) {
        cerr << usage << endl;
        exit(1);
    }

    const int n{atoi(argv[1])};
    const int nEvectors{atoi(argv[2])};
    const string sortByString{argv[3]};

    SortByType sortBy{};
    if (sortByString == "real") {
        sortBy = SortByType::real;
    } else if (sortByString == "abs") {
        sortBy = SortByType::abs;
    } else {
        cerr << usage << endl;
        exit(2);
    }

    // Class with the function that returns A * x for given x.
    VectorFunction& vf(VectorFunction::instance(n));

    // Make a random initial vector.
    VectorXd initial{n};
    for (int i = {0}; i < n; ++i) {
        initial(i) = randDouble(-0.5, 0.5);
    }

    cout << "Initial vector: (" << initial.transpose() << ")" << endl;

    eig_ns::ArnoldiEig<VectorXd> eig{
        std::bind(&VectorFunction::f, &vf, _1),
        vectorIP,
        initial,
        n,
        sortBy,
        true
    };

    // Compute the correct eigenvalues and eigenvectors.
    const MatrixXd A{vf.getMatrix()};
    Eigen::EigenSolver<MatrixXd> eigenSolver{A, true};

    VectorXcd correct_evalues{eigenSolver.eigenvalues()};
    MatrixXcd correct_evectors{eigenSolver.eigenvectors()};
    eig.sortEigPair(correct_evalues, correct_evectors);

    cout << "\nCorrect eigenvalues =" << endl;
    cout << correct_evalues << endl;
    cout << "\nCorrect eigenvectors =" << endl;
    cout << correct_evectors << endl;

    // Compute only eigenvalues.
    cout << "\nChecking eigenvalues alone." << endl;
    VectorXcd eValues = eig.eVals();
    cout << "\nEigenvalues =" << endl;
    cout << eValues << endl;

    // Check eigenvalues and eigenvectors together.
    cout << "\nChecking eigenvalues and eigenvectors together." << endl;
    vector<vector<VectorXd>> v(
        n,
        vector<VectorXd>(2, VectorXd::Zero(n))
    );
    eig.eig(eValues, v, nEvectors);
    // Convert from vector<VectorXd> to MatrixXcd.
    MatrixXcd eVectors{n, nEvectors};
    for (int i = {0}; i < n; ++i) {
        for (int j = {0}; j < nEvectors; ++j) {
            eVectors(i, j) = std::complex<double>{v[j][0](i), v[j][1](i)};
        }
    }

    cout << "\nEigenvalues =" << endl;
    cout << eValues << endl;
    cout << "\nEigenvectors = " << endl;
    cout << eVectors << endl;

    // Check that the Arnoldi eigenvectors are correct.
    cout << "\nArnoldi eigenvectors divided by correct eigenvectors ="
         << endl;
    cout << eVectors.array() / correct_evectors.array().leftCols(nEvectors)
         << endl;
}
