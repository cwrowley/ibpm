/// \file AnalyzeEig.H
/// \brief A template for evaluating and printing eigenvalues and eigenvectors.
/// \author Kevin K. Chen (Princeton University)

#ifndef ANALYZE_EIG_H
#define ANALYZE_EIG_H

#include <functional>
#include <vector>
#include <Eigen>
#include <Core>


#include "templateAliases.h"

// ************************************************************************** //

namespace eig_ns {

using std::function;
using std::vector;

using Eigen::VectorXcd;

// ************************************************************************** //

/// Indicator for whether the dynamical operator is continuous or discrete.
enum class OperatorType {continuous, discrete};

// ************************************************************************** //
// AnalyzeEig template declaration.

/// \brief Class for evaluating and printing eigenvalues and eigenvectors along
/// equilibria arcs, e.g., as required by Continuation.

/// \tparam Vec Vector type for \c AnalyzeEig.
template <typename Vec>
class AnalyzeEig {
private:
    /// \brief \c std::function that takes in a state, an equilibrium point, and
    /// a parameter, and returns a new state.
    using LinearFunctionLam =
        function<Vec (const Vec&, const Vec&, double)>;
    /// \brief \c std::function that takes in a \c Function and computes its
    /// eigenvalues and eigenvectors.
    using EigSolver =
        function<void (alias::Function<Vec>, VectorXcd&, vector<vector<Vec>>&)>;
    /// \c std::function that prints eigenvalues.
    using EValuesPrintFunction = function<void (const VectorXcd&, double)>;
    /// \c std::function that prints eigenvectors.
    using EVectorsPrintFunction =
        function<void (const vector<vector<Vec>>&, double)>;

public:
    /// Constructor.
    AnalyzeEig(
        const EigSolver& eig_solver,
        const LinearFunctionLam& eig_stepper, // Function for eigensolver.
        OperatorType eig_stepper_type,
        const EValuesPrintFunction evalues_print,
        const EVectorsPrintFunction evectors_print,
        const Vec& x_eq,
        double& lam
    );
    ~AnalyzeEig() noexcept {} ///< Destructor (empty).

    /// Compute the eigenvalues and eigenvectors at the given equilibrium.
    int analyzeEig(const Vec& x);

private:
    const EigSolver eig_solver_; ///< Eigenvalue solver.
    const LinearFunctionLam eig_stepper_; ///< Function for the eigensolver.
    const OperatorType eig_stepper_type_; ///< \c continuous or \c discrete.
    const EValuesPrintFunction evalues_print_; ///< Eigenvalues printer.
    const EVectorsPrintFunction evectors_print_; ///< Eigenvectors printer.

    Vec x_eq_; ///< \brief The current root.
    double& lam_; /// Reference to the continuation parameter.

    AnalyzeEig(const AnalyzeEig&) = delete;
    AnalyzeEig(AnalyzeEig&&) = delete;

    void operator=(const AnalyzeEig&) = delete;
    void operator=(AnalyzeEig&&) = delete;

    /// Compute \f$\mathrm{eig_stepper\_} (\mathbf{x}; \mathbf{x}_e,
    /// \lambda)\f$.
    Vec eigStepperAtLam(const Vec& x) const;

    /// Return the number of entries that have non-negative real part.
    static int numPos(const VectorXcd& x);
    /// Return the number of entries with unity or greater magnitude.
    static int numAbsGT1(const VectorXcd& x);
};

// ************************************************************************** //
// AnalyzeEig template member functions.

/// \param[in] eig_solver An \c std::function that computes eigenvalues and
/// eigenvectors from a function.
/// \param[in] eig_stepper An \c std::function representing the operator to
/// eigendecompose.
/// \param[in] eig_stepper_type Either \c OperatorType::continuous or \c
/// OperatorType::discrete, describing the way eig_stepper handles time.
/// \param[in] evalues_print An \c std::function that prints eigenvalues.
/// \param[in] evectors_print An \c std::function that prints eigenvectors.
/// \param[in] x_eq Any \c Vec object, for initialization.
/// \param[in] lam A reference to the continuation parameter.
template <typename Vec>
AnalyzeEig<Vec>::AnalyzeEig(
    const EigSolver& eig_solver,
    const LinearFunctionLam& eig_stepper,
    OperatorType eig_stepper_type,
    const EValuesPrintFunction evalues_print,
    const EVectorsPrintFunction evectors_print,
    const Vec& x_eq,
    double& lam
) :
    eig_solver_{eig_solver},
    eig_stepper_{eig_stepper},
    eig_stepper_type_{eig_stepper_type},
    evalues_print_{evalues_print},
    evectors_print_{evectors_print},

    x_eq_{x_eq},
    lam_(lam)
{
}

/// After the eigenvalues and eigenvectors have been solved, print them.
///
/// \param[in] x An equilibrium state.
/// \return The number of unstable eigenvalues.
template <typename Vec>
int AnalyzeEig<Vec>::analyzeEig(const Vec& x) {
    x_eq_ = x;

    VectorXcd eValues{};
    vector<vector<Vec>> eVectors;

    using namespace std::placeholders;
    // Compute the eigendecomposition.
    eig_solver_(
        std::bind(&AnalyzeEig::eigStepperAtLam, this, _1),
        eValues,
        eVectors
    );

    evalues_print_(eValues, lam_);
    evectors_print_(eVectors, lam_);

    if (eig_stepper_type_ == OperatorType::continuous) {
        return numPos(eValues);
    } else {
        return numAbsGT1(eValues);
    }
}

/// This form is necessary, since the standard function interface has exactly
/// one input vector and one output vector.  This form "hides" the value of the
/// parameter from the interface.
///
/// \param[in] x Input vector.
/// \return fx Output vector, equal to \f$\mathrm{eig_stepper\_}
/// (\mathbf{x}; \mathbf{x}_e, \lambda)\f$.
template <typename Vec>
inline Vec AnalyzeEig<Vec>::eigStepperAtLam(const Vec& x) const {
    return eig_stepper_(x, x_eq_, lam_);
}

/// \param[in] x A complex vector.
/// \return The number of entries \f$x_j\f$ in \f$\mathbf{x}\f$ with
/// \f$\mathrm{Re} (x_j) \ge 0\f$.
/// \sa numAbsGT1
template <typename Vec>
int AnalyzeEig<Vec>::numPos(const VectorXcd& x) {
    int n{0};

    for (int i = {0}; i < x.size(); ++i) {
        if (x[i].real() >= 0) {
            ++n;
        }
    }

    return n;
}

/// \param[in] x A complex vector.
/// \return The number of entries \f$x_j\f$ in \f$\mathbf{x}\f$ with \f$|x_j|
/// \ge 1\f$.
/// \sa numPos
template <typename Vec>
int AnalyzeEig<Vec>::numAbsGT1(const VectorXcd& x) {
    int n = 0;

    for (int i = 0; i < x.size(); ++i) {
        if (std::abs(x[i]) >= 1.) {
            ++n;
        }
    }

    return n;
}

// ************************************************************************** //

} // End namespace eig_ns.

#endif // ANALYZE_EIG_H
