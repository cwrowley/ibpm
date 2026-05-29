/// \file LinearSolver.h
/// \brief Interface for solving linear systems \f$\mathbf{A x} = \mathbf{b}\f$.
/// \author Kevin K. Chen, Clarence W. Rowley (Princeton University)

#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <functional>

#include "templateAliases.h"

// ************************************************************************** //

/// The namespace used for all linear-solver-related constructs.
namespace linear_solver_ns {

// ************************************************************************** //
// LinearSolver template declaration.

/// \brief Template interface for solving linear systems \f$\mathbf{A x} =
/// \mathbf{b}\f$.

/// This interface is intended for large problems where the matrix
/// \f$\mathbf{A}\f$ is not explicitly known.  Instead, the linear system
/// solving is assumed to require only a function that computes the product
/// \f$\mathbf{A x}\f$ from the vector \f$\mathbf{x}\f$.
///
/// \tparam Vec Vector type for \c LinearSolver.
/// \sa Gmres
template <typename Vec>
class LinearSolver {
public:
    LinearSolver(bool verbose = false); ///< Constructor.
    virtual ~LinearSolver() noexcept {} ///< Destructor (empty).

    /// Solve \f$\mathbf{A x} = \mathbf{b}\f$ for \f$\mathbf{x}\f$.
    virtual Vec solve(const alias::Function<Vec>& a, const Vec& b) const = 0;

protected:
    const bool verbose_; ///< \c true to activate verbose printing.

private:
    LinearSolver(const LinearSolver&) = delete;
    LinearSolver(LinearSolver&&) = delete;

    void operator=(const LinearSolver&) = delete;
    void operator=(LinearSolver&&) = delete;
};

// ************************************************************************** //
// LinearSolver template member functions.

template <typename Vec>
LinearSolver<Vec>::LinearSolver(bool verbose) : verbose_{verbose} {
}

// ************************************************************************** //

} // End namespace linear_solver_ns

#endif // LINEARSOLVER_H
