/// \file NewtonArmijo.h
/// \brief A template for using Newton's method with the Armijo line search
/// method, for nonlinear root-finding.
/// \author Kevin K. Chen, Clarence W. Rowley, and Zhanhua Ma (Princeton
/// University)

#ifndef NEWTON_ARMIJO_H
#define NEWTON_ARMIJO_H

#include <iostream>

#include "LinearSolver.h"
#include "Newton.h"
#include "templateAliases.h"

// ************************************************************************** //

namespace newton_ns {

// ************************************************************************** //
// NewtonArmijo template declaration.

/// \brief A template for using Newton's method with the Armijo line search
/// method, for nonlinear root-finding.

/// This template is a child class of the \c Newton template, and is essentially
/// identical except for one key difference: instead of using the update rule
/// \f[
///     \mathbf{x}_{j+1} = \mathbf{x}_j - \mathbf{h},
/// \f]
/// where \f$\mathbf{h}\f$ is computed using a linear systems solver, it uses
/// \f[
///     \mathbf{x}_{j+1} = \mathbf{x}_j - \lambda \mathbf{h},
/// \f]
/// where \f$\lambda\f$ is determined by the Armijo rule.
///
/// The default values of \f$\sigma_0\f$, \f$\sigma_1\f$,
/// \f$\mathrm{max_arm_iter}\f$, and \f$\alpha\f$ are the same as C. T. Kelley's
/// nsoli.m.
///
/// Reference: C. T. Kelley, \a "Solving Nonlinear Equations with Newton's
/// Method," and Kelley's MATLAB code nsoli.m.
///
/// \sa Newton
template <typename Vec>
class NewtonArmijo : public Newton<Vec> {
public:
    /// Constructor when \f$\mathbf{f}\f$ is known at construction time.
    NewtonArmijo(
        const alias::Function<Vec>& f,
        const LinearSolver<Vec>& linearSolver,
        const alias::Norm<Vec>& norm,
        double tol = Newton<Vec>::tol_default_,
        int max_iter = Newton<Vec>::max_iter_default_,
        double jacobian_dx = Newton<Vec>::jacobian_dx_default_,
        int jacobian_order = Newton<Vec>::jacobian_order_default_,
        int max_arm_iter = max_arm_iter_default_,
        double s0 = s0Default_,
        double s1 = s1Default_,
        double a = aDefault_,
        bool verbose = false
    );
    /// Constructor when \f$\mathbf{f}\f$ is not known at construction time.
    NewtonArmijo(
        const LinearSolver<Vec>& linearSolver,
        const alias::Norm<Vec>& norm,
        double tol = Newton<Vec>::tol_default_,
        int max_iter = Newton<Vec>::max_iter_default_,
        double jacobian_dx = Newton<Vec>::jacobian_dx_default_,
        int jacobian_order = Newton<Vec>::jacobian_order_default_,
        int max_arm_iter = max_arm_iter_default_,
        double s0 = s0Default_,
        double s1 = s1Default_,
        double a = aDefault_,
        bool verbose = false
    );

    ~NewtonArmijo() noexcept override {}

private:
    const int max_arm_iter_; ///< The maximum number of Armijo iterations.
    const double s0_; ///< The value of \f$\sigma_0\f$.
    const double s1_; ///< The value of \f$\sigma_1\f$.
    const double a_; ///< The value of \f$\alpha\f$.

    /// Default \c max_arm_iter_.
    static constexpr int max_arm_iter_default_{10};
    static constexpr double s0Default_{0.1}; ///< Default \f$\sigma_0\f$.
    static constexpr double s1Default_{0.5}; ///< Default \f$\sigma_1\f$.
    static constexpr double aDefault_{1.e-4}; ///< Default \f$\alpha\f$.

    NewtonArmijo(const NewtonArmijo&) = delete;
    NewtonArmijo(NewtonArmijo&&) = delete;

    void operator=(const NewtonArmijo&) = delete;
    void operator=(NewtonArmijo&&) = delete;

    /// Return the step multiplier \f$\lambda\f$ computed from the Armijo rule.
    double stepMultiplier(
        const alias::Function<Vec>& f,
        const Vec& x,
        const Vec& fx,
        const Vec& h
    ) const override;

    /// Apply a 3-point safeguarded parabolic line search minimization.
    double parab3p(
        double lam_curr,
        double lam_prev,
        double fx_norm2,
        double f_curr_norm2,
        double f_prev_norm2
    ) const;
    /// Get the threshold for quitting the Armijo line search.
    double threshold(double lambda, double fx_norm) const;

    /// Print convergence information.
    void printConvergence(
        int iter,
        double lambda,
        double f_curr_norm,
        double tol
    ) const;
    void checkInputs() const; ///< Check for valid inputs.
};

// ************************************************************************** //
// NewtonArmijo template member functions.

/// \param[in] f An \c std::function that computes an output vector from an
/// input vector and a parameter.
/// \param[in] linearSolver An \f$\mathbf{A x} = \mathbf{b}\f$ solver.
/// \param[in] norm An \c std::function that returns the norm of a \c Vec.
/// \param[in] tol The norm tolerance below which Newton quits.
/// \param[in] max_iter The maximum number of iterations before quitting.
/// \param[in] jacobian_dx The step size for Jacobian-vector products.
/// \param[in] jacobian_order The order of accuracy for Jacobian-vector
/// products.
/// \param[in] max_arm_iter The maximum number of Armijo iterations to run.
/// \param[in] s0 The value of \f$\sigma_0\f$.
/// \param[in] s1 The value of \f$\sigma_1\f$.
/// \param[in] a The value of \f$\alpha\f$.
/// \param[in] verbose \c true to active verbose printing.
template <typename Vec>
NewtonArmijo<Vec>::NewtonArmijo(
    const alias::Function<Vec>& f,
    const LinearSolver<Vec>& linearSolver,
    const alias::Norm<Vec>& norm,
    double tol,
    int max_iter,
    double jacobian_dx,
    int jacobian_order,
    int max_arm_iter,
    double s0,
    double s1,
    double a,
    bool verbose
) :
    Newton<Vec>{
        f,
        linearSolver,
        norm,
        tol,
        max_iter,
        jacobian_dx,
        jacobian_order,
        verbose
    },
    max_arm_iter_{max_arm_iter},
    s0_{s0},
    s1_{s1},
    a_{a}
{
    checkInputs();
}

/// \param[in] linearSolver An \f$\mathbf{A x} = \mathbf{b}\f$ solver.
/// \param[in] norm An \c std::function that returns the norm of a \c Vec.
/// \param[in] tol The norm tolerance below which Newton quits.
/// \param[in] max_iter The maximum number of iterations before quitting.
/// \param[in] jacobian_dx The step size for Jacobian-vector products.
/// \param[in] jacobian_order The order of accuracy for Jacobian-vector
/// products.
/// \param[in] max_arm_iter The maximum number of Armijo iterations to run.
/// \param[in] s0 The value of \f$\sigma_0\f$.
/// \param[in] s1 The value of \f$\sigma_1\f$.
/// \param[in] a The value of \f$\alpha\f$.
template <typename Vec>
NewtonArmijo<Vec>::NewtonArmijo(
    const LinearSolver<Vec>& linearSolver,
    const alias::Norm<Vec>& norm,
    double tol,
    int max_iter,
    double jacobian_dx,
    int jacobian_order,
    int max_arm_iter,
    double s0,
    double s1,
    double a,
    bool verbose
) :
    NewtonArmijo<Vec>{
        std::bind(&Newton<Vec>::unset, _1),
        linearSolver,
        norm,
        tol,
        max_iter,
        jacobian_dx,
        jacobian_order,
        max_arm_iter,
        s0,
        s1,
        a,
        verbose
    }
{
}

/// \param[in] f An \c std::function that computes an output vector from an
/// input vector and a parameter.
/// \param[in] x The state variable at the start of the step.
/// \param[in] fx The value of \f$\mathbf{f} (\mathbf{x})\f$.
/// \param[in] h The step to take by Newton's method.
/// \return The scalar factor \f$\lambda\f$ by which to modify the Newton step
/// vector.
template <typename Vec>
double NewtonArmijo<Vec>::stepMultiplier(
    const alias::Function<Vec>& f,
    const Vec& x,
    const Vec& fx,
    const Vec& h
) const {
    using std::cout;
    using std::cerr;
    using std::endl;

    Vec xCurr{x};
    double lambda{1.};
    xCurr -= h * lambda;

    Vec fCurr{f(xCurr)};

    const double fx_norm{this->norm_(fx)};
    const double fx_norm2{fx_norm * fx_norm};
    double f_curr_norm{this->norm_(fCurr)};
    double f_curr_norm2{f_curr_norm * f_curr_norm};
    double f_prev_norm2{f_curr_norm2};

    double lam_curr{lambda};
    double lam_prev{1.};

    double tol{threshold(lambda, fx_norm)};

    if (this->verbose_) {
        printConvergence(0, lambda, f_curr_norm, tol);
    }

    for (int iter = {0}; iter < max_arm_iter_ && f_curr_norm >= tol; ++iter) {
        if (this->verbose_) {
            cout << "\nArmijo iteration " << iter + 1 << "." << endl;
        }

        if (iter == 0) {
            lambda *= s1_;
        } else { // Apply the three-point parabolic model.
            lambda = parab3p(
                lam_curr,
                lam_prev,
                fx_norm2,
                f_curr_norm2,
                f_prev_norm2
            );
        }

        // Update x.
        xCurr = x - h * lambda;
        lam_prev = lam_curr;
        lam_curr = lambda;

        // Update function evaluation and norms.
        fCurr = f(xCurr);
        f_prev_norm2 = f_curr_norm2;
        f_curr_norm = this->norm_(fCurr);
        f_curr_norm2 = f_curr_norm * f_curr_norm;

        // Update the threshold.
        tol = threshold(lambda, fx_norm);

        if (this->verbose_) {
            printConvergence(iter + 1, lambda, f_curr_norm, tol);
        }
    }

    if (f_curr_norm < threshold(lambda, fx_norm)) {
        if (this->verbose_) {
            cout << "\n||f(x + lambda * h)|| is below threshold.  "
                 << "Stopping Armijo line search.\n" << endl;
        }
    } else {
        cerr << "Armijo line search failed to converge after " << max_arm_iter_
             << " iterations.\n" << endl;
    }

    return lambda;
}

/// This estimates the value of \f$\lambda\f$ leading to the minimum of
/// \f$\|\mathbf{f} (\mathbf{x} + \lambda \mathbf{h})\|^2\f$.  It is safeguarded
/// in the sense that the algorithm must return \f$\lambda \in [\sigma_0
/// \lambda_{\mathrm{curr}}, \sigma_1 \lambda_{\mathrm{curr}}]\f$.
///
/// \param[in] lam_curr The current step size \f$\lambda_{\mathrm{curr}}\f$.
/// \param[in] lam_prev the previous step size \f$\lambda_{\mathrm{prev}}\f$.
/// \param[in] fx_norm2 \f$\|\mathbf{f} (\mathbf{x})\|^2\f$.
/// \param[in] f_curr_norm2 \f$\|\mathbf{f} (\mathbf{x} +
/// \lambda_{\mathrm{curr}} \mathbf{h})\|^2\f$.
/// \param[in] f_prev_norm2 \f$\|\mathbf{f} (\mathbf{x} +
/// \lambda_{\mathrm{prev}} \mathbf{h})\|^2\f$.
/// \return The minimizing value of \f$\lambda\f$ given by the parabolic model.
template <typename Vec>
double NewtonArmijo<Vec>::parab3p(
    double lam_curr,
    double lam_prev,
    double fx_norm2,
    double f_curr_norm2,
    double f_prev_norm2
) const {
    double lambda;
    double c2{
        lam_prev * (f_curr_norm2 - fx_norm2) -
        lam_curr * (f_prev_norm2 - fx_norm2)
    };

    if (c2 >= 0) {
        lambda = s1_ * lam_curr;
    } else {
        double c1{
            lam_curr * lam_curr * (f_prev_norm2 - fx_norm2) -
            lam_prev * lam_prev * (f_curr_norm2 - fx_norm2)
        };
        lambda = -c1 * 0.5 / c2;

        if (lambda < s0_ * lam_curr) {
            lambda = s0_ * lam_curr;
        } else if (lambda > s1_ * lam_curr) {
            lambda = s1_ * lam_curr;
        }
    }

    return lambda;
}

/// \param[in] lambda The current value of \f$\lambda\f$.
/// \param[in] fx_norm The value of \f$\|\mathbf{f} (\mathbf{x})\|\f$.
/// \return \f$(1 - \alpha \lambda) \|\mathbf{f} (\mathbf{x})\|\f$.
template <typename Vec>
inline double NewtonArmijo<Vec>::threshold(
    double lambda,
    double fx_norm
) const {
    return (1 - a_ * lambda) * fx_norm;
}

/// \param[in] iter The iteration number.
/// \param[in] lambda The current value of \f$\lambda\f$.
/// \param[in] f_curr_norm The value of \f$\|\mathbf{f} (\mathbf{x} + \lambda
/// \mathbf{h})\|\f$.
/// \param[in] tol The tolerance on \c f_curr_norm for quitting the iterator.
template <typename Vec>
inline void NewtonArmijo<Vec>::printConvergence(
    int iter,
    double lambda,
    double f_curr_norm,
    double tol
) const {
    std::cout << "\nArmijo iteration " << iter << ":" << " lambda = " << lambda
              << "; ||f(x + lambda * h)|| = " << f_curr_norm << "; tolerance = "
              << tol << "." << std::endl;
}

/// Throw an \c InvalidValue exception if certain inputs don't make sense.
template <typename Vec>
void NewtonArmijo<Vec>::checkInputs() const {
    using aux_ns::InvalidValue;

    if (max_arm_iter_ < 0) {
        throw InvalidValue{"max_arm_iter >= 0."};
    }
    if (s0_ >= s1_) {
        throw InvalidValue{"s0 < s1 required."};
    }
    if (a_ <= 0.) {
        throw InvalidValue{"a > 0 required."};
    }
    if (a_ >= 1.) {
        throw InvalidValue{"a < 1 required."};
    }
}

// ************************************************************************** //

} // End namespace newton_ns

#endif // NEWTON_ARMIJO_H
