/// \file Continuation.h
/// \brief A template for performing numerical continuation on a smooth arc of
/// equilibria.
///
/// Includes the \c Continuation template and \c Param struct.
/// \author Kevin K. Chen (Princeton University)

#ifndef CONTINUATION_H
#define CONTINUATION_H

#include <complex>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "Newton.h"
#include "exceptions.h"

// ************************************************************************** //

/// \namespace continuation_ns
/// \brief The namespace used for all numerical-continuation-related constructs.
namespace continuation_ns {

using std::cout;
using std::endl;
using std::function;
using std::string;

using newton_ns::Newton;

// ************************************************************************** //
// Continuation template declaration.

/// \brief A template for performing continuation on a smooth arc of equilibria,
/// using extrapolation.

/// Given
///   - a function \f$\mathbf{f} (\mathbf{x}, \lambda)\f$,
///   - a root solver (currently restricted to \c newton_ns::Newton<Vec>),
///   - an eigenvalue solver,
///   - and a number of miscellaneous parameters,
///
/// this template solves for the solution arc \f$\mathbf{x}_0 (\lambda)\f$ such
/// that \f$\mathbf{f} (\mathbf{x}_0 (\lambda), \lambda) = \mathbf{0}\f$.  For
/// each value of \f$\lambda\f$, it performs continuation by using extrapolation
/// as a predictor and the root solver as the corrector.  The first
/// extrapolation is 0th order, and all subsequent extrapolations are 1st order.
///
/// After an equilibrium has been found with the root solver, the eigenvalue
/// solver is run to report the found eigenvalues and in particular, then number
/// of unstable eigenvalues.  If the number of unstable eigenvalues has changed
/// since the last root found, a message will announce this.
///
/// This template uses adaptive parameter stepping.  It will always first try to
/// take the largest parameter step.  Every time the Newton solver throws a \c
/// NewtonError, \c Continuation will halve the parameter step size.  If this
/// results in a step size below the given minimum step size, then \c
/// Continuation gives up and throws a \c ContinuationError.
///
/// \tparam Vec Vector type for \c Continuation.
template <typename Vec>
class Continuation {
private:
    /// \brief \c std::function that takes in a state and a parameter, and
    /// computes a new state.
    using FunctionLam = function<Vec (const Vec&, double)>;
    /// \c std::function that prints a state with a \c double parameter.
    using PrintFunction = function<void (const Vec&, double)>;
    /// \brief \c std::function that analyzes the eigenvalues at an input state
    /// and returns the number of unstable eigenvalues.
    using AnalyzeEigFunction = function<int (const Vec&)>;

public:
    class Param; ///< Class for holding continuation parameter details.

    /// Constructor.
    Continuation(
        const FunctionLam& f,
        Newton<Vec>& newton,
        const Vec& x,
        const PrintFunction& print,
        const AnalyzeEigFunction& analyze_eig,
        const AnalyzeEigFunction& analyze_adjoint_eig,
        Param& param
    );
    ~Continuation() noexcept {} ///< Destructor (empty).

    void solve(); ///< Compute the equilibrium arc.

private:
    const FunctionLam f_; ///< Function we want to solve the solution arc of.
    Newton<Vec>& newton_; ///< Reference to Newton solver object.

    Vec x_; ///< Current estimated root location.
    Vec xPrev_; ///< Previous root location.
    const PrintFunction print_; ///< Function to print each solution found.

    /// Perform eigendecomposition.
    const AnalyzeEigFunction analyze_eig_;
    /// Perform adjoint eigendecomposition.
    const AnalyzeEigFunction analyze_adjoint_eig_;

    Param& param_; ///< Reference to continuation parameters object.
    int n_unstable_eig_prev_{-1}; ///< Previous number of unstable eigenvalues.
    int extrapolate_order_{0}; ///< 0 on first iteration; 1 afterwards.

    Continuation(const Continuation&) = delete;
    Continuation(Continuation&&) = delete;

    void operator=(const Continuation&) = delete;
    void operator=(Continuation&&) = delete;

    /// Compute \f$\mathbf{f} (\mathbf{x}, \lambda)\f$.
    Vec fAtLam(const Vec& x) const;
    bool finished() const; /// \brief Check if continuation is finished.

    void step(); ///< Extrapolate to the next \f$\lambda\f$, and solve Newton.
    void extrapolate(); ///< Extrapolate to the next \f$\lambda\f$.
    void solveRoot(); ///< Solve the Newton iterator.
    void computeEig(); ///< Compute direct and adjoint eigendecomposition.
    void dlamReduce(); ///< Reduce step size when the Newton iterator fails.

    static double sign(double x); ///< Signum function.
};

// ************************************************************************** //
// Param class declaration.

/// \brief Data container for holding values related to the continuation
/// parameter \f$\lambda\f$.

/// \tparam Vec Vector type for \c Continuation.
/// \sa Continuation<Vec>::Continuation
template <typename Vec>
class Continuation<Vec>::Param {
public:
    const string name; ///< Name of continuation parameter (e.g., \c "lambda").
    double lam; ///< Current value of \f$\lambda\f$
    const double lam_f; ///< Final value of \f$\lambda\f$
    const double dlam_max; ///< Maximum (in magnitude) \f$\lambda\f$ step.
    const double dlam_min; ///< Minimum (in magnitude) \f$\lambda\f$ step.
    const double dlam_sign; ///< Sign of \c dlam_max and \c dlam_min.
    double dlam; ///< Current Continuation parameter step \f$\Delta \lambda\f$.
    double dlam_prev; ///< Previous parameter step \f$\Delta \lambda\f$.

    /// Constructor.
    Param(
        const string& name,
        double lam,
        double lam_f,
        double dlam_max,
        double dlam_min
    );
    ~Param() noexcept {} ///< Destructor (empty).

private:
    Param(const Param&) = delete;
    Param(Param&&) = delete;

    void operator=(const Param&) = delete;
    void operator=(Param&&) = delete;

    void checkInputs() const; ///< Check for valid inputs.
};

// ************************************************************************** //
// Continuation template member functions.

/// \param[in] f An \c std::function that computes an output vector from an
/// input vector and a parameter.
/// \param[in] newton A Newton root-solver object.
/// \param[in] x Initial guess for the equilibrium at \f$\lambda\f$.
/// \param[in] print An \c std::function that prints the state and the
/// parameter, e.g., when a root is found.
/// \param[in] analyze_eig An \c std::function that computes the
/// eigendecomposition at an input vector, and returns the number of unstable
/// eigenvalues.
/// \param[in] analyze_adjoint_eig An \c std::function that computes the adjoint
/// eigendecomposition at an input vector, and returns the number of unstable
/// eigenvalues.
/// \param[in] param A \c Param struct containing continuation parameters.
///
/// \sa Param, AnalyzeEig<Vec>::analyzeEig(const Vec&)
template <typename Vec>
Continuation<Vec>::Continuation(
    const FunctionLam& f,
    Newton<Vec>& newton,
    const Vec& x,
    const PrintFunction& print,
    const AnalyzeEigFunction& analyze_eig,
    const AnalyzeEigFunction& analyze_adjoint_eig,
    Param& param
) :
    f_{f},
    newton_(newton),

    x_{x},
    xPrev_{x},
    print_{print},

    analyze_eig_{analyze_eig},
    analyze_adjoint_eig_{analyze_adjoint_eig},

    param_(param)
{
    using namespace std::placeholders;
    // The Newton object cannot take a function interface with a lam parameter.
    // Give it a function where lam is used without being explicitly passed in.
    newton_.set_f(std::bind(&Continuation<Vec>::fAtLam, this, _1));
}

// Compute the equilibrium arc.
template <typename Vec>
void Continuation<Vec>::solve() {
    // First, solve for the root at the initial value of lam, starting at the
    // given initial condition.
    try {
        solveRoot();
    } catch (newton_ns::NewtonError& ne) {
        std::ostringstream s{};
        s << "\nNewton failed: " << ne.what() << "\n"
          << "Error: Newton failed to converge during the initial Newton "
          << "solve.  Exiting continuation.";
        throw ContinuationError{s.str()};
    }

    computeEig();

    // Then, extrapolate to the next lam and solve for the root there.
    while (!finished()) {
        step();
    }

    cout << "\n*****" << endl;
    cout << "\nContinuation completed." << endl;
}

/// This form is necessary, since the standard function interface has exactly
/// one input vector and one output vector.  This form "hides" the value of the
/// parameter from the interface.
///
/// \param[in] x Input vector.
/// \return fx Output vector, equal to \f$\mathbf{f} (\mathbf{x},
/// \lambda)\f$.
template <typename Vec>
inline Vec Continuation<Vec>::fAtLam(const Vec& x) const {
    return f_(x, param_.lam);
}

/// The continuation is finished if \f$\lambda\f$ is at or past the final value.
///
/// \return Boolean value.
template <typename Vec>
inline bool Continuation<Vec>::finished() const {
    return param_.dlam_sign * (param_.lam_f - param_.lam) <= 0.;
}

/// If the Newton iteration fails, then reduce the magnitude of the step size
/// \f$\Delta \lambda\f$ and try again.
template <typename Vec>
void Continuation<Vec>::step() {
    // Remember these values, in case we need to reset to them.
    Vec x{x_};
    Vec xPrev{xPrev_};

    extrapolate();

    try {
        solveRoot();
    } catch (newton_ns::NewtonError& ne) {
        cout << "\nNewton failed: " << ne.what() << endl;

        // Reset these values.
        x_ = x;
        xPrev_ = xPrev;

        dlamReduce(); // Try a smaller lam.
        return;
    }

    // If we reached this point, then Newton succeeded.  Reset dlam.
    param_.dlam_prev = param_.dlam;
    param_.dlam = param_.dlam_max;

    computeEig();

    // Set extrapolation order to 1 after the first iteration.
    extrapolate_order_ = 1;
}

/// Use 0th order (constant) extrapolation the first time, but switch to 1st
/// order (linear) extrapolation for all subsequent calls.
template <typename Vec>
void Continuation<Vec>::extrapolate() {
    cout << "\n**********" << endl;
    cout << "\nExtrapolating to " << param_.name << " = "
         << param_.lam + param_.dlam << "." << endl;
    param_.lam += param_.dlam;

    xPrev_ = x_;

    // Note: for 0th order extrapolation, we don't need to do anything more.
    if (extrapolate_order_ == 1) { // Linear extrapolation.
        x_ += (x_ - xPrev_) * (param_.dlam / param_.dlam_prev);
    } else if (extrapolate_order_ != 0) {
        throw aux_ns::OrderError{
            "Extrapolation orders above 1 are not implemented."
        };
    }
}

/// Also, print the state and parameter value after the root has been found.
template <typename Vec>
void Continuation<Vec>::solveRoot() {
    cout << "\nBegin Newton at " << param_.name << " = " << param_.lam << "."
         << endl;
    newton_.solve(x_); // Ideally returns 0.
    print_(x_, param_.lam);
}

/// Print a message if a bifurcation has been found.
template <typename Vec>
void Continuation<Vec>::computeEig() {
    cout << "\nComputing the eigendecomposition of the linearization." << endl;
    int nUnstableEig{analyze_eig_(x_)};

    // On the first call to solveRoot(), n_unstable_eig_prev_ = -1.  Otherwise,
    // nUnstableEig >= 0.
    if (n_unstable_eig_prev_ >= 0 && nUnstableEig != n_unstable_eig_prev_) {
        cout << "\n***** Bifurcation detected." << endl;
        cout << n_unstable_eig_prev_ << " unstable eigenvalues at "
             << param_.name << " = " << param_.lam - param_.dlam << "." << endl;
        cout << nUnstableEig << " unstable eigenvalues at " << param_.name
             << " = " << param_.lam << "." << endl;
    }
    n_unstable_eig_prev_ = nUnstableEig;

    cout << "\nComputing the adjoint eigendecomposition of the linearization."
         << endl;
    analyze_adjoint_eig_(x_);
}

/// Try reducing the parameter step size by a factor of 2.  If this results in a
/// step size below the allowable range, then give up and throw a \c
/// ContinuationError exception.
template <typename Vec>
void Continuation<Vec>::dlamReduce() {
    using std::abs;

    // Go back to the last working value of param_.lam.
    param_.lam -= param_.dlam;
    param_.dlam *= 0.5;

    // If we already reduced dlam as much as we could, then give up.
    if (abs(param_.dlam) < abs(param_.dlam_min)) {
        std::ostringstream s{};
        s << "Error: Newton failed to converge for the allowable range of "
          << param_.name << ".  Exiting continuation.";
        throw ContinuationError{s.str()};
    }

    cout << "Trying " << param_.name << " = " << param_.lam + param_.dlam << "."
         << endl;
}

/// \param[in] x A real scalar.
/// \return
/// \f[
///     \mathrm{sgn} (x) =
///     \left\{
///         \begin{array}{rl}
///         -1, & x < 0 \\
///         0, & x = 0 \\
///         1, & x > 0
///         \end{array}
///     \right.
/// \f]
template <typename Vec>
inline double Continuation<Vec>::sign(double x) {
    return x > 0. ? 1. : (x < 0. ? -1. : 0.);
}

// ************************************************************************** //
// Param class member functions.

/// \param[in] name_in The name of the parameter, e.g., \c "lambda".
/// \param[in] lam_in The initial value of \f$\lambda\f$, i.e., at the beginning
/// of the solution arc.
/// \param[in] lam_f_in The final value of \f$\lambda\f$, i.e., at the end of
/// the solution arc.
/// \param[in] dlam_max_in The largest value of \f$|\Delta \lambda|\f$ to use
/// when advancing to the next value of \f$\lambda\f$.
/// \param[in] dlam_min_in The smallest value of \f$|\Delta \lambda|\f$ to use
/// when advancing to the next value of \f$\lambda\f$.
template <typename Vec>
Continuation<Vec>::Param::Param(
    const string& name_in,
    double lam_in,
    double lam_f_in,
    double dlam_max_in,
    double dlam_min_in
) :
    name{name_in},
    lam{lam_in},
    lam_f{lam_f_in},
    dlam_max{dlam_max_in},
    dlam_min{dlam_min_in},
    dlam_sign{sign(dlam_max)},
    dlam{dlam_max}
{
    checkInputs();
}

/// Throw an \c InvalidValue exception if certain inputs don't make sense.
template <typename Vec>
void Continuation<Vec>::Param::checkInputs() const {
    using std::abs;
    using std::ostringstream;
    using aux_ns::InvalidValue;

    if (dlam_min == 0) {
        ostringstream s{};
        s << "dlam_min = " << dlam_min << " must be nonzero.";
        throw InvalidValue{s.str()};
    }
    if (abs(dlam_max) < abs(dlam_min)) {
        ostringstream s{};
        s << "abs(dlam_max) = " << abs(dlam_max)
          << " cannot be smaller than abs(dlam_min) = " << abs(dlam_min)
          << ".";
        throw InvalidValue{s.str()};
    }
    if (sign(dlam_max) != sign(dlam_min)) {
        throw InvalidValue{"dlam_max and dlam_min must have the same sign."};
    }
}

// ************************************************************************** //

} // End namespace continuation_ns

#endif // CONTINUATION_H
