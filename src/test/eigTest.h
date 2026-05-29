#ifndef EIG_TEST_H
#define EIG_TEST_H

#include <Eigen>
#include <Core>
#include "random.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// ************************************************************************** //
// VectorFunction class declaration.

class VectorFunction {
public:
    static VectorFunction& instance(int n); // Get the singleton instance.

    const MatrixXd& getMatrix() const {return a_;}
    // Function evaluation.
    VectorXd f(const VectorXd& x) const {return a_ * x;}

private:
    const int n_; // The dimension of the vector function input and output.
    MatrixXd a_; // The linear transformation from input to output.

    explicit VectorFunction(int n);
    VectorFunction(const VectorFunction&) = delete;
    VectorFunction(VectorFunction&&) = delete;

    ~VectorFunction() noexcept {}

    void operator=(const VectorFunction&) = delete;
    void operator=(VectorFunction&&) = delete;
};

// ************************************************************************** //
// VectorFunction member functions.

// Get the singleton instance.
inline VectorFunction& VectorFunction::instance(int n) {
    static VectorFunction vf{n};
    return vf;
}

// Constructor.
inline VectorFunction::VectorFunction(int n) : n_(n), a_(n, n) {
    for (int i = {0}; i < n_; ++i) {
        for (int j = {0}; j < n_; ++j) {
            a_(i, j) = randDouble(-5., 5.);
        }
    }
}

// ************************************************************************** //
// Other functions.

// VectorXd inner product interface.
inline double vectorIP(const VectorXd& x, const VectorXd& y) {
    return y.dot(x);
}

// Return the input.  This is needed for EmpiricalEig's conversion to/from
// VectorXd and Vec, the latter of which is just VectorXd in this case.
inline VectorXd identity(const VectorXd& x) {
    return x;
}

// ************************************************************************** //

#endif // EIG_TEST_H
