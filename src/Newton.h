/// \file Newton.h
/// \brief A template for using Newton's method to perform nonlinear
/// \modified for periodic orbit computing with unkown period line 263 to line 281
/// root-finding.
/// \author Kevin K. Chen, Clarence W. Rowley (Princeton University)

#ifndef NEWTON_H
#define NEWTON_H

#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <iostream>

#include "Jacobian.h"
#include "LinearSolver.h"
#include "exceptions.h"
#include "templateAliases.h"

// ************************************************************************** //

/// The namespace used for all Newton-related constructs.
namespace newton_ns {

using namespace std::placeholders;

using alias::Function;
using alias::Norm;
using linear_solver_ns::LinearSolver;

using std::ostream;
using std::ofstream;

// ************************************************************************** //
 //Newton template declaration.

 //A template for using Newton's method to perform nonlinear root-finding.

 //Newton's method attempts to solve \f$\mathbf{f} (\mathbf{x}) = \mathbf{0}\f$
 //by iterating with
 //\f[
 //    \mathbf{x}_{j+1} =
 //    \mathbf{x}_j -
 //    \left(
 //        \left.
 //            \frac{\partial \mathbf{f}}{\partial \mathbf{x}}
 //        \right|_{\mathbf{x}_j}
 //    \right)^{-1}
 //    \mathbf{f} (\mathbf{x}_j).
 //\f]
 //If the vector \f$\mathbf{x}\f$ is large, then inverting the Jacobian
 //\f$\left.\partial \mathbf{f} / \partial \mathbf{x}\right|_{\mathbf{x}_j}\f$
 //is not practical, so instead we use the two-step process
/*
//\f{eqnarray*}{
 //    \left.
 //        \frac{\partial \mathbf{f}}{\partial \mathbf{x}}
 //    \right|_{\mathbf{x}_j}
 //    \mathbf{h} &=&
 //    \mathbf{f} (\mathbf{x}_j) \\
 //    \mathbf{x}_{j+1} &=& \mathbf{x}_j - \mathbf{h}.
 //\f}
 */
 //We solve the first equation using an iterative solver, often based on a
 //Krylov method.  Several different alternatives are possible for this
 //iterative solver (e.g. GMRES, Bi-CGSTAB, TFQMR; see Kelley, "Iterative
 //Methods for Linear and Nonlinear Equations"), so here we decouple the
 //details of the linear solver from the Newton iteration.

 //This design uses \c std::function objects both for the nonlinear function
 //\f$\mathbf{f}\f$, and for the Jacobian \f$\left.\partial \mathbf{f} /
 //\partial \mathbf{x}\right|_{\mathbf{x}_j}\f$ (which here is a function
 //object).

 //\tparam Vec Vector type for \c Newton.
 //\sa NewtonArmijo
template <typename Vec>
class Newton {
private:
    using FuncAtIter = std::function<void (const Vec&, int)>;

public:
    /// Constructor when \f$\mathbf{f}\f$ is known at construction time.
    Newton(
        const Function<Vec>& f,
        const LinearSolver<Vec>& linear_solver,
        const Norm<Vec>& norm,
        double tol = tol_default_,
        int max_iter = max_iter_default_,
        double jacobian_dx = jacobian_dx_default_,
        int jacobian_order = jacobian_order_default_,
        bool verbose = false
    );
    /// Constructor when \f$\mathbf{f}\f$ is not known at construction time.
    Newton(
        const LinearSolver<Vec>& linear_solver,
        const Norm<Vec>& norm,
        double tol = tol_default_,
        int max_iter = max_iter_default_,
        double jacobian_dx = jacobian_dx_default_,
        int jacobian_order = jacobian_order_default_,
        bool verbose = false
    );

    virtual ~Newton() noexcept {} ///< Destructor (empty).

    /// Solve Newton's method and return \f$\mathbf{f}(\mathbf{x})\f$.
    Vec solve(Vec& x) const;

    /// Set the function to solve the root of.
    void set_f(const Function<Vec>& f);
    /// Set the function to be run at the end of every iteration.
    void setFuncAtIter(const FuncAtIter& f);

protected:
    /// Function that returns the norm of the input \c Vec.
    const Norm<Vec> norm_;
    const bool verbose_; ///< \c true to activate verbose printing.

    static constexpr int max_iter_default_{10}; ///< Default maximum iterations.
    static constexpr double tol_default_{1.e-6}; ///< Default tolerance.

    /// Default Jacobian order.
    static constexpr int jacobian_order_default_{1};
    static constexpr double jacobian_dx_default_{1.e-3}; ///< Default step.

    static Vec unset(const Vec& x); ///< An unset \c Function.

private:
    /// Function that returns \f$\mathbf{f} (\mathbf{x})\f$.
    Function<Vec> f_;
    const LinearSolver<Vec>& linear_solver_; ///< \brief Linear systems solver.
    /// Function to be run after every iteration.
    FuncAtIter func_at_iter_{[](const Vec&, int) {}};

    const double tol_; ///< Norm tolerance below which Newton quits.
    const int max_iter_; ///< The maximum number of iterations before quitting.
    const double jacobian_dx_; ///< The step size for Jacobian-vector products.
    const int jacobian_order_; ///< Accuracy order for Jacobian-vector products.

    Newton(const Newton&) = delete;
    Newton(Newton&&) = delete;

    void operator=(const Newton&) = delete;
    void operator=(Newton&&) = delete;

    void step(Vec& x, const Vec& fx) const; ///< \brief Take one Newton step.
    /// The factor to modify the Newton step size with.
    virtual double stepMultiplier(
        const Function<Vec>& f,
        const Vec& x,
        const Vec& fx,
        const Vec& h
    ) const;
    /// Print convergence information.
    void printConvergence(int iter, double norm_fx) const;
    void checkInputs() const; ///< Check for valid inputs.
};

// ************************************************************************** //
// Newton template member functions.

/// \param[in] f An \c std::function that computes an output vector from an
/// input vector and a parameter.
/// \param[in] linear_solver An \f$\mathbf{A x} = \mathbf{b}\f$ solver.
/// \param[in] norm An \c std::function that returns the norm of a \c Vec.
/// \param[in] tol The norm tolerance below which the Newton iterator quits.
/// \param[in] max_iter The maximum number of iterations before quitting.
/// \param[in] jacobian_dx The step size for Jacobian-vector products.
/// \param[in] jacobian_order The order of accuracy for Jacobian-vector
/// products.
/// \param[in] verbose \c true to active verbose printing.
template <typename Vec>
Newton<Vec>::Newton(
    const Function<Vec>& f,
    const LinearSolver<Vec>& linear_solver,
    const Norm<Vec>& norm,
    double tol,
    int max_iter,
    double jacobian_dx,
    int jacobian_order,
    bool verbose
) :
    norm_(norm),
    verbose_{verbose},
    f_(f),
    linear_solver_(linear_solver),
    tol_(tol),
    max_iter_(max_iter),
    jacobian_dx_(jacobian_dx),
    jacobian_order_(jacobian_order)
{
    checkInputs();
}

/// This constructor is useful, for instance, when constructing a \c
/// Continuation object, where \f$\mathbf{f}\f$ depends on the continuation
/// parameter.
///
/// \param[in] linear_solver An \f$\mathbf{A x} = \mathbf{b}\f$ solver.
/// \param[in] norm An \c std::function that returns the norm of a \c Vec.
/// \param[in] tol The norm tolerance below which the Newton iterator quits.
/// \param[in] max_iter The maximum number of iterations before quitting.
/// \param[in] jacobian_dx The step size for Jacobian-vector products.
/// \param[in] jacobian_order The order of accuracy for Jacobian-vector
/// products.
/// \param[in] verbose \c true to active verbose printing.
template <typename Vec>
Newton<Vec>::Newton(
    const LinearSolver<Vec>& linear_solver,
    const Norm<Vec>& norm,
    double tol,
    int max_iter,
    double jacobian_dx,
    int jacobian_order,
    bool verbose
) :
    Newton<Vec>{
        std::bind(&Newton<Vec>::unset, _1),
        linear_solver,
        norm,
        tol,
        max_iter,
        jacobian_dx,
        jacobian_order,
        verbose
    }
{
}

/// The solver is run until \f$\|\mathbf{f} (\mathbf{x})\| < \mathrm{tol}\f$, or
/// the maximum number of iterations is reached.  If the maximum number of
/// iterations is reached without convergence, then a \c NewtonError is thrown.
/// A \c NewtonError is also thrown if <tt>step(Vec& x, const Vec& fx)</tt>
/// throws an exception, e.g., if the linear solver fails.
///
/// \param[in,out] x The initial guess of the solution; modified in place to
/// return the final estimate as well.
/// \return The value of \f$\mathbf{f} (\mathbf{x})\f$, which should ideally be
/// close to \f$\mathbf{0}\f$.
template <typename Vec>
Vec Newton<Vec>::solve(Vec& x) const {
    using std::cout;
    using std::endl;
    using std::ostringstream;

    Vec fx{f_(x)};
    double norm_fx{norm_(fx)};

    if (verbose_) {
        printConvergence(0, norm_fx);
    }

    int iter{};
    for (iter = 0; iter < max_iter_ && norm_fx > tol_; ++iter) {
        if (verbose_) {
            cout << "\nNewton iteration " << iter + 1 << "." << endl;
        }

		// compute phase condition (period) only used for periodic orbit with unknown
        // period
		int PhaseCal_flag = 0;
		// save PhaseCal_flag
		ofstream file;
		file.open("../examples/PhaseCal_flag.txt");
		file << PhaseCal_flag;
		file.close();
		cout << "PhaseCal_flag: " << PhaseCal_flag << endl;
		//
		f_(x); // Compute f(x) to get phase condition
		PhaseCal_flag = PhaseCal_flag + 1;
		// save updated PhaseCal_flag
		ofstream file_PF;
		file_PF.open("../examples/PhaseCal_flag.txt");
		file_PF << PhaseCal_flag;
		file_PF.close();
		cout << "PhaseCal_flag: " << PhaseCal_flag << endl;
		//

        // Take a Newton step, updating x.
        try {
            step(x, fx);
        } catch (linear_solver_ns::LinearSolverError& lse) {
            ostringstream s{};
            s << "Linear solver failed: " << lse.what();
            throw NewtonError{s.str()};
        } catch (std::exception& e) {
            ostringstream s{};
            s << "Newton step failed: " << e.what();
            throw NewtonError{s.str()};
        }

        fx = f_(x); // Compute f(x).
        norm_fx = norm_(fx);

        // Extra function to be run after each iteration, if the user specified
        // one.
        func_at_iter_(x, iter + 1);

        if (verbose_) {
            printConvergence(iter + 1, norm_fx);
        }
    }

    if (norm_fx <= tol_) {
        if (verbose_) {
            cout << "\n||fx|| is at or below tolerance.  "
                 << "Stopping Newton iteration.\n" << endl;
        }
    } else {
        ostringstream s{};
        s << "Newton failed to converge after " << iter << " iterations.  "
          << "||fx|| = " << norm_fx << " is above tolerance " << tol_
          << ".";
        throw NewtonError{s.str()};
    }

    return fx;
}

/// \param[in] f The std::function to solve the root of.
template <typename Vec>
inline void Newton<Vec>::set_f(const Function<Vec>& f) {
    f_ = f;
}

/// \param[in] f The std::function that should be run after every iteration.
/// Note that the second input to f is a zero-index iteration number.
template <typename Vec>
inline void Newton<Vec>::setFuncAtIter(const FuncAtIter& f) {
    func_at_iter_ = f;
}

/// \param[in] x The input vector.
/// \return fx The output vector.
template <typename Vec>
inline Vec Newton<Vec>::unset(const Vec& x) {
    throw(std::logic_error{"Newton::f_ is not set."});

    // Should not reach this point; put this here to turn off compiler warnings.
    return x;
}

/// The update algorithm is
/// \f[
///     \mathbf{x} \leftarrow
///     \mathbf{x} -
///     \left(
///         \left.
///             \frac{\partial \mathbf{f}}{\partial \mathbf{x}}
///         \right|_{\mathbf{x}}
///     \right)^{-1}
///     \mathbf{f} (\mathbf{x}).
/// \f]
///
/// \param[in,out] x The state at the start of the Newton step; returned as the
/// new state the stepper iterated to.
/// \param[in] fx The value of \f$\mathbf{f} (\mathbf{x})\f$.
template <typename Vec>
inline void Newton<Vec>::step(Vec& x, const Vec& fx) const {
    using std::cout;
    using std::endl;

    const Jacobian<Vec> Df_x{f_, x, fx, jacobian_dx_, jacobian_order_};

    if (verbose_) {
        cout << "\nSolving Df(x) h = f(x)." << endl;
    }
    Vec h{linear_solver_.solve(Df_x, fx)};

    if (verbose_) {
        cout << "\nTaking Newton step." << endl;
    }
    // The stepMultiplier is just 1 for the Newton method; it can be a different
    // scalar for modified Newton methods.
    x -= h * stepMultiplier(f_, x, fx, h);
}

/// For the Newton method, this value is always 1.  Child classes may return
/// other values.
///
/// \param[in] f An std::function that computes an output vector from an input
/// vector and a parameter.  Unused.
/// \param[in] x The state variable at the start of the step.  Unused.
/// \param[in] fx The value of \f$\mathbf{f} (\mathbf{x})\f$.  Unused.
/// \param[in] h The step to take by Newton's method.  Unused.
/// \return 1.0.
template <typename Vec>
inline double Newton<Vec>::stepMultiplier(
    const Function<Vec>&,
    const Vec&,
    const Vec&,
    const Vec&
) const {
    return 1.;
}

/// \param[in] iter The iteration number.
/// \param[in] norm_fx The value of \f$||\mathbf{f} (\mathbf{x})||\f$.
template <typename Vec>
inline void Newton<Vec>::printConvergence(int iter, double norm_fx) const {
    std::cout << "\nNewton iteration " << iter << ": ||fx|| = " << norm_fx
              << "; tolerance = " << tol_ << "." << std::endl;
}

/// Throw an \c InvalidValue exception if certain inputs don't make sense.
template <typename Vec>
void Newton<Vec>::checkInputs() const {
    using aux_ns::InvalidValue;

    if (tol_ < 0) {
        throw InvalidValue{"tol >= 0 required."};
    }
    if (max_iter_ < 0) {
        throw InvalidValue{"max_iter >= 0 required."};
    }
    if (jacobian_dx_ <= 0) {
        throw InvalidValue{"jacobian_dx > 0 required."};
    }
}

// ************************************************************************** //

} // End namespace newton_ns

#endif // NEWTON_H
