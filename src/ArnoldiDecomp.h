/// \file ArnoldiDecomp.h
/// \brief A template for performing Arnoldi decomposition.
/// \author Kevin K. Chen, Jonathan H. Tu (Princeton University)

#ifndef ARNOLDI_DECOMP_H
#define ARNOLDI_DECOMP_H

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

#include <Eigen>
#include <Core>

#include "exceptions.h"
#include "templateAliases.h"

// ************************************************************************** //

/// \namespace arnoldi_ns
/// \brief The namespace used for all Arnoldi-related constructs.
namespace arnoldi_ns {

using std::function;
using std::vector;

using Eigen::MatrixXd;
using Eigen::VectorXcd;

// ************************************************************************** //
// ArnoldiDecomp template declaration.

/// A template for performing Arnoldi decomposition.

/// We approximate the action of an operator A using the decomposition
/// \f[
///     \mathbf{A Q} = \mathbf{Q}^\mathrm{T} \mathbf{h}.
/// \f]
/// Recall that the eigenvalues of \f$\mathbf{H}(1:n, 1:n)\f$ (i.e.,
/// \f$\mathbf{H}\f$ minus its last row) are approximations to the eigenvalues
/// of \f$\mathbf{A}\f$.  We say that \f$\mathbf{H}\f$ is "full" when its values
/// in the \f$n\f$th column have been filled in, ie when we can compute \f$n\f$
/// eigenvalues.
///
/// \tparam Vec Vector type for \c ArnoldiDecomp.
/// \sa ArnoldiEig
template <typename Vec>
class ArnoldiDecomp {
public:
    /// Constructor.
    ArnoldiDecomp(
        const alias::Function<Vec>& a,
        const alias::InnerProduct<Vec>& ip,
        const Vec& b,
        int n,
        bool verbose = false
    );
    ~ArnoldiDecomp() noexcept {} ///< Destructor (empty).

    void step(); ///< Take one iteration of the Arnoldi decomposition.
    void step(int n_steps); ///< Compute \c n_steps more iterations.
    void stepTo(int n_steps); ///< Iterate until \c currentSize() is n_steps.
    void stepAll(); ///< Iterate until \c currentSize() = \c maxSize().

    /// Return complex linear combinations of the columns of \f$\mathbf{Q}\f$.
    vector<Vec> qTimes(const VectorXcd& y);
    /// Return real linear combinations of the columns of \f$\mathbf{Q}\f$.
    Vec qTimes(const vector<double>& y);

    int currentSize() const; ///< Return the size of the Hessenberg matrix.
    int maxSize() const; ///< Return the maximum Arnoldi size.
    MatrixXd getH() const; ///< Return the Hessenberg matrix.
    /// Return the latest assigned column of Hessenberg matrix.
    MatrixXd getHLastCol() const;
    vector<Vec> getQ() const; ///< Return the matrix of Arnoldi vectors.

private:
    /// Function that returns \f$\mathbf{A x}\f$.
    const alias::Function<Vec> a_;
    const alias::InnerProduct<Vec> ip_; ///< Inner product function.

    const int n_;///< Maximum number of approximating eigenvalues.
    int iter_; ///< Iteration number.

    MatrixXd h_; ///< Rectangular Hessenberg matrix.
    vector<Vec> q_{}; ///< Matrix of Arnoldi vectors.
    Vec v_; ///< Arnoldi vector.
    const double tol_{1.e-10}; ///< Tolerance below which an exception is raised.

    const bool verbose_; ///< \c true to activate verbose printing.

    ArnoldiDecomp(const ArnoldiDecomp&) = delete;
    ArnoldiDecomp(ArnoldiDecomp&&) = delete;

    void operator=(const ArnoldiDecomp&) = delete;
    void operator=(ArnoldiDecomp&&) = delete;
};

// ************************************************************************** //
// ArnoldiDecomp template member functions.

/// \param[in] a The \c std::function for the products of \f$\mathbf{A}\f$.
/// \param[in] ip The \c std::function defining the inner product for \c Vec.
/// \param[in] b The initial guess for the Arnoldi decomposition.
/// \param[in] n The dimension of the approximation to \f$\mathbf{A}\f$.
/// \param[in] verbose \c true to activate verbose printing.
template <typename Vec>
ArnoldiDecomp<Vec>::ArnoldiDecomp(
    const alias::Function<Vec>& a,
    const alias::InnerProduct<Vec>& ip,
    const Vec& b,
    int n,
    bool verbose
) :
    a_{a},
    ip_{ip},

    n_{n},
    iter_{0},

    h_{MatrixXd::Zero(n + 1, n)},
    v_{b},

    verbose_{verbose}
{
    q_.push_back(b / sqrt(ip_(b, b)));
}

/// Compute one more element of \f$\mathbf{Q}\f$, and the corresponding
/// row and column of \f$\mathbf{H}\f$.
///
/// \sa step(int n_steps), stepTo(int n_steps), and stepAll()
template <typename Vec>
void ArnoldiDecomp<Vec>::step() {
    if (verbose_) {
        ++iter_;
        std::cout << "\nArnoldi iteration " << iter_ << "." << std::endl;
    }

    // Only continue to step if h_ isn't full.
    int k{currentSize()};
    if (k < n_) {
        v_ = a_(q_[k]);

        for (int j = {0}; j <= k; ++j) {
            h_(j, k) = ip_(v_, q_[j]);
            v_ -= q_[j] * h_(j, k);
        }

        h_(k + 1, k) = sqrt(ip_(v_, v_));

        // Check the last element at the lastest assigned column, unless the
        // lastest column is the last column of h_.
        if (k < n_ - 1 && h_(k + 1, k) < tol_) {
            std::ostringstream s{};
            s << "h_(k + 1, k) = " << h_(k + 1, k) << " is below tolerance.";
            throw ArnoldiError{s.str()};
        }

        if (h_(k + 1, k) == 0) {
            q_.push_back(
                v_ * std::numeric_limits<double>::quiet_NaN()
            );
        } else {
            q_.push_back(v_ / h_(k + 1, k));
        }

        if (verbose_) {
            std::cout << "\nCurrent Hessenberg size is " << currentSize() << "."
                      << std::endl;
        }
    }
}

/// Compute \c n_steps more elements of \f$\mathbf{Q}\f$, until \c currentSize()
/// is \f$n\f$.
///
/// \param[in] n_steps Number of iterations to take.
/// \sa step(), stepTo(int n_steps), and stepAll()
template <typename Vec>
inline void ArnoldiDecomp<Vec>::step(int n_steps) {
    for (int n = {0}; n < n_steps; ++n) {
        step();
    }
}

/// \param[in] n_steps Arnoldi decomposition size to iterate to.
/// \sa step(), step(int n_steps), and stepAll()
template <typename Vec>
inline void ArnoldiDecomp<Vec>::stepTo(int n_steps) {
    step(n_steps - currentSize());
}

/// \sa step(), step(int n_steps), and stepTo(int n_steps)
template <typename Vec>
inline void ArnoldiDecomp<Vec>::stepAll() {
    step(n_ - currentSize());
}

/// This first method treats the general case of a complex vector
/// \f$\mathbf{y}\f$.
///
/// \param[in] y A complex vector.
/// \return The product \f$\mathbf{Q y}\f$, given as an \c std::vector where the
/// first element is the real part and the second element is the imaginary part.
/// \sa qTimes(const vector<double>& y)
template <typename Vec>
inline vector<Vec> ArnoldiDecomp<Vec>::qTimes(const VectorXcd& y) {
    assert(y.size() == static_cast<int>(q_.size()) - 1);
    vector<Vec> x(2, q_[0] * 0.);

    for (int i = {0}; i < y.size(); ++i) {
        x[0] += q_[i] * y(i).real();
        x[1] += q_[i] * y(i).imag();
    }

    return x;
}

/// This second method treats the simple case of a real vector \f$\mathbf{y}\f$.
///
/// \param[in] y A real vector.
/// \return The product \f$\mathbf{Q y}\f$.
/// \sa qTimes(const VectorXcd& y)
template <typename Vec>
inline Vec ArnoldiDecomp<Vec>::qTimes(const vector<double>& y) {
    assert(y.size() == q_.size() - 1);
    Vec x{q_[0] * 0.};

    for (unsigned i = {0}; i < y.size(); ++i) {
        x += q_[i] * y[i];
    }

    return x;
}

/// Note that due to the algorithm setup, there will be one more vector in
/// \f$\mathbf{Q}\f$ than there will be eigenvalues for \f$\mathbf{H}_n\f$; it
/// starts from 0.
///
/// \return The size of the Hessenberg matrix \f$\mathbf{H}_n\f$.
/// \sa maxSize()
template <typename Vec>
inline int ArnoldiDecomp<Vec>::currentSize() const {
    return static_cast<int>(q_.size()) - 1;
}

/// \return The maximum number of approximating eigenvalues.
/// \sa currentSize()
template <typename Vec>
inline int ArnoldiDecomp<Vec>::maxSize() const {
    return n_;
}

/// \return The Hessenberg matrix \f$\mathbf{H}_n\f$.
/// \sa getHLastCol()
template <typename Vec>
inline MatrixXd ArnoldiDecomp<Vec>::getH() const {
    return h_;
}

/// \return The latest assigned column of the Hessenberg matrix
/// \f$\mathbf{H}_n\f$.
/// \sa getH()
template <typename Vec>
inline MatrixXd ArnoldiDecomp<Vec>::getHLastCol() const {
    return h_.col(currentSize() - 1);
}

/// \return The matrix \f$\mathbf{Q}\f$ of Arnoldi vectors.
/// \sa getH()
template <typename Vec>
inline vector<Vec> ArnoldiDecomp<Vec>::getQ() const {
    return q_;
}

//************************************************************************** //

} // End namespace arnoldi_ns

#endif // ARNOLDI_DECOMP_H
