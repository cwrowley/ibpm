/// \file Jacobian.h
/// \brief A template for computing Jacobian products.
/// \author Kevin K. Chen, Clarence W. Rowley (Princeton University)

#ifndef JACOBIAN_H
#define JACOBIAN_H

#include <functional>

#include "exceptions.h"
#include "templateAliases.h"

// ************************************************************************** //

namespace newton_ns {

// ************************************************************************** //
// Jacobian template declaration.

/// A template for computing the Jacobian-vector products.

/// The product this template computes is
/// \f[
///     \mathbf{y} = \left.
///         \frac{\partial \mathbf{f}}{\partial \mathbf{x}}
///     \right|_{\mathbf{x}_0}
///     \mathbf{x}.
/// \f]
/// The computation is done by finite differencing.  The algorithm for
/// first-order finite differencing is
/// \f[
///     \left.
///         \frac{\partial \mathbf{f}}{\partial \mathbf{x}}
///     \right|_{\mathbf{x}_0}
///     \mathbf{x}
///     =
///     \frac{
///         \mathbf{f} (\mathbf{x}_0 + \epsilon \mathbf{x}) -
///         \mathbf{f} (\mathbf{x}_0)
///     }
///     {\epsilon} +
///     \mathcal{O} (\epsilon)
/// \f]
/// and the \f$\epsilon\f$-order error is truncated.  This computation requires
/// one function evaluation if \f$\mathbf{f} (\mathbf{x}_0)\f$ is given to the
/// constructor, and two evaluations otherwise.
///
/// The algorithm for second-order finite differencing is
/// \f[
///     \left.
///         \frac{\partial \mathbf{f}}{\partial \mathbf{x}}
///     \right|_{\mathbf{x}_0}
///     \mathbf{x}
///     =
///     \frac{
///         \mathbf{f} (\mathbf{x}_0 + \epsilon \mathbf{x}) -
///         \mathbf{f} (\mathbf{x}_0 - \epsilon \mathbf{x})
///     }
///     {2 \epsilon} +
///     \mathcal{O} (\epsilon^2)
/// \f]
/// and the \f$\epsilon^2\f$-order error is truncated.  This computation
/// requires two function evaluations.
///
/// \tparam Vec Vector type for \c Jacobian.
template <typename Vec>
class Jacobian {
public:
    /// Constructor if \f$\mathbf{f} (\mathbf{x}_0)\f$ is known.
    Jacobian(
        const alias::Function<Vec>& f,
        const Vec& x0,
        const Vec& fx0,
        double eps,
        int order = order_default_
    );
    /// Constructor if \f$\mathbf{f} (\mathbf{x}_0)\f$ is not known.
    Jacobian(
        const alias::Function<Vec>& f,
        const Vec& x0,
        double eps,
        int order = order_default_
    );
    Jacobian(const Jacobian&) = default; ///< Copy constructor.
    Jacobian(Jacobian&&) = default; ///< Move constructor.

    ~Jacobian() noexcept {} ///< Destructor (empty).

    /// Compute \f$\mathbf{y} = \mathbf{J x}\f$ using finite differencing.
    Vec operator()(const Vec& x) const;

private:
    /// Function that returns \f$\mathbf{f} (\mathbf{x})\f$.
    const alias::Function<Vec> f_;
    const Vec x0_; ///< State where the Jacobian is evaluated.
    Vec fx0_;  ///< Value of \f$\mathbf{f} (\mathbf{x}_0)\f$.

    const double eps_; ///< Step size for finite differencing.
    const int order_; ///< Order of accuracy for finite-differencing (1 or 2).
    static constexpr int order_default_{1}; ///< Default \c order_.

    void checkInputs() const; ///< Check for valid inputs.

    void operator=(const Jacobian&) = delete;
    void operator=(Jacobian&&) = delete;
};

// ************************************************************************** //
// Jacobian template member functions.

/// \param[in] f An \c std::function that this template computes the Jacobian
/// of.
/// \param[in] x0 The equilibrium point about which the Jacobian is computed.
/// \param[in] fx0 The value of \f$\mathbf{f} (\mathbf{x}_0)\f$, ideally but
/// possibly not exactly \f$\mathbf{0}\f$.
/// \param[in] eps The step size used for finite differncing.
/// \param[in] order The order of accuracy for finite-differencing (1 or 2
/// only).  Note that \c fx0 is unused if \c order is 2.
/// \sa Jacobian(const alias::Function<Vec>& f, const Vec& x0, double eps)
template <typename Vec>
Jacobian<Vec>::Jacobian(
    const alias::Function<Vec>& f,
    const Vec& x0,
    const Vec& fx0,
    double eps,
    int order
) :
    f_{f},
    x0_{x0},
    fx0_{fx0},
    eps_{eps},
    order_{order}
{
    checkInputs();
}

/// \param[in] f An \c std::function that this template computes the Jacobian
/// of.
/// \param[in] x0 The equilibrium point about which the Jacobian is computed.
/// \param[in] eps The step size used for finite differncing.
/// \param[in] order The order of accuracy for finite-differencing (1 or 2
/// only).
/// \sa Jacobian(const alias::Function<Vec>& f, const Vec& x0, const Vec& fx0,
/// double eps)
template <typename Vec>
Jacobian<Vec>::Jacobian(
    const alias::Function<Vec>& f,
    const Vec& x0,
    double eps,
    int order
) :
    Jacobian{f, x0, f(x0), eps, order}
{
}

/// \param[in] x The vector to multiply the Jacobian by.
/// \return The product \f$\mathbf{J x}\f$.
template <typename Vec>
Vec Jacobian<Vec>::operator()(const Vec& x) const {
    Vec y{f_(x0_ + x * eps_)};

    if (order_ == 1) {
        y -= fx0_;
        y /= eps_;
    } else { // order_ == 2
        Vec temp{f_(x0_ - x * eps_)};
        y -= temp;
        y /= 2 * eps_;
    }

    return y;
}

/// Throw an \c InvalidValue exception if certain inputs don't make sense.
template <typename Vec>
void Jacobian<Vec>::checkInputs() const {
    using aux_ns::InvalidValue;

    if (eps_ == 0.) {
        throw InvalidValue{"eps != 0 required."};
    }
    if (order_ != 1 && order_ != 2) {
        throw InvalidValue{"order must be 1 or 2."};
    }
}

// ************************************************************************** //

} // End namespace newton_ns

#endif // JACOBIAN_H
