/// \file Gmres.h
/// \brief A template for performing the Generalized Minimal Residual method.
/// \author Kevin K. Chen, Clarence W. Rowley (Princeton University)

#ifndef GMRES_H
#define GMRES_H

#include <cassert>
#include <cmath>

#include <iostream>
#include <vector>

#include <Eigen>
#include <Core>

#include "ArnoldiDecomp.h"
#include "LinearSolver.h"
#include "exceptions.h"
#include "templateAliases.h"

// ************************************************************************** //

namespace linear_solver_ns {

using alias::Function;
using alias::InnerProduct;
using std::vector;
using Eigen::MatrixXd;

// ************************************************************************** //
// Gmres template declaration.

/// A template for performing the Generalized Minimal Residual method.

/// GMRES is one of a number of iterative algorithms for estimating the solution
/// \f$\mathbf{x}\f$ to the linear system \f$\mathbf{A x} = \mathbf{b}\f$.  It
/// is intended for systems where the operator \f$\mathbf{A}\f$ is so large that
/// its explicit form is generally unavailable.  Instead, it only requires the
/// function that returns \f$\mathbf{A x}\f$ given \f$\mathbf{x}\f$.
///
/// \tparam Vec Vector type for \c Gmres.
/// \sa LinearSolver
template <typename Vec>
class Gmres : public LinearSolver<Vec> {
public:
    /// Constructor.
    Gmres(
        const InnerProduct<Vec>& ip,
        int maxIter,
        double tol = 1.e-6,
        bool verbose = false
    );
    ~Gmres() noexcept override {} ///< Destructor (empty).

    /// Solve \f$\mathbf{A x} = \mathbf{b}\f$ for \f$\mathbf{x}\f$.
    Vec solve(const Function<Vec>& a, const Vec& b) const override;

    class GivensRotation; // Class for Givens rotation products.

    const InnerProduct<Vec> ip_; ///< Inner product function.
    const int maxIter_; ///< \brief Maximum number of iterations to perform.
    /// \brief Tolerance on \f$\|\mathbf{A x} - \mathbf{b}\| / \|\mathbf{b}\|\f$
    /// for quitting.
    const double tol_;

    Gmres(const Gmres&) = delete;
    Gmres(Gmres&&) = delete;

    void operator=(const Gmres&) = delete;
    void operator=(Gmres&&) = delete;

    /// Upper-triangulize the currently-assigned part of \f$\mathbf{H}\f$.
    void upperTriangulize(
        vector<GivensRotation>& givensList,
        MatrixXd& H,
        int Hsize
    ) const;
    /// \brief Solve \f$\mathbf{H y} = \mathbf{g}\f$ by back substitution,
    /// assuming \f$\mathbf{H}\f$ is upper triangular.
    vector<double> back_substitute(
        const MatrixXd& H,
        const vector<double>& g,
        int Hsize
    ) const;
    /// Print convergence information.
    void printConvergence(int iter, double rho, double errtol) const;
    void checkInputs() const; ///< Check for valid inputs.
};

// ************************************************************************** //
// GivensRotation class declaration.

/// Class for Givens rotation products.

/// \tparam Vec Vector type for \c GivensRotation.
template <typename Vec>
class Gmres<Vec>::GivensRotation {
public:
    /// Construct a Givens Rotation matrix.
    GivensRotation(int k, double c, double s);
    GivensRotation(GivensRotation&&) = default;

    ~GivensRotation() noexcept {} ///< Destructor (empty).

    /// Apply to the \f$i\f$th column of a matrix \f$\mathbf{H}\f$.
    void applyToACol(MatrixXd& H, int i);
    /// Apply to a column vector \f$\mathbf{v}\f$.
    void apply(vector<double>& v);

private:
    int k_; ///< Index (with 1-indexing) where the 2 x 2 block is.
    double c_; ///< Upper left and lower right element of the 2 x 2 block.
    double s_; ///< Lower left element of the 2 x 2 block.

    GivensRotation(const GivensRotation&) = delete;

    void operator=(const GivensRotation&) = delete;
    void operator=(GivensRotation&&) = delete;

    /// Apply a two-dimensional rotation.
    void twoDRotation(double& v1, double& v2);
};

// ************************************************************************** //
// Gmres template member functions.

/// \param[in] ip The \c std::function defining the inner product for \c Vec.
/// \param[in] maxIter The maximum number of iterations to run GMRES.
/// \param[in] tol The value of \f$\|\mathbf{A x} - \mathbf{b}\| /
/// \|\mathbf{b}\|\f$ below which GMRES should quit.
template <typename Vec>
Gmres<Vec>::Gmres(
    const InnerProduct<Vec>& ip,
    int maxIter,
    double tol,
    bool verbose
) :
    LinearSolver<Vec>{verbose},
    ip_{ip},
    maxIter_{maxIter},
    tol_{tol}
{
    checkInputs();
}

/// \param[in] A The \c std::function for products of \f$\mathbf{A}\f$.
/// \param[in] b The right-hand side of the linear systems equation.
/// \return x The estimated solution.
template <typename Vec>
Vec Gmres<Vec>::solve(const Function<Vec>& a, const Vec& b) const
{
    using std::cout;
    using std::cerr;
    using std::endl;

    double rho{std::sqrt(ip_(b, b))};
    double errtol{tol_ * rho};

    arnoldi_ns::ArnoldiDecomp<Vec> arnoldi{a, ip_, b, maxIter_};
    MatrixXd H{MatrixXd::Zero(maxIter_ + 1, maxIter_)};
    vector<double> g(1, rho); // Intermediate vector; initialized as rho.
    vector<GivensRotation> givensList; // Givens rotation series.

    if (this->verbose_) {
        printConvergence(0, rho, errtol);
    }

    for (int iter{0}; iter < maxIter_ && rho > errtol; ++iter) {
        if (this->verbose_) {
            cout << "\nGMRES iteration " << iter + 1 << "." << endl;
        }
        try {
            arnoldi.step();
        } catch (arnoldi_ns::ArnoldiError& ae) {
            std::ostringstream s{};
            s << "Arnoldi iteration failed: " << ae.what();
            throw LinearSolverError{s.str()};
        }

        // Assign the iter-th column.
        H.block(0, iter, maxIter_ + 1, 1) = arnoldi.getHLastCol();
        // Upper-triangularize H.
        upperTriangulize(givensList, H, arnoldi.currentSize());

        g.push_back(0);
        assert(static_cast<int>(g.size()) == iter + 2);
        givensList[iter].apply(g);

        rho = std::abs(g[iter+1]);

        if (this->verbose_) {
            printConvergence(iter + 1, rho, errtol);
        }
    }

    if (this->verbose_) {
        if (rho <= errtol) {
            cout << "\nrho is at or below tolerance " << errtol
                 << ".  Stopping GMRES iteration.\n" << endl;
        }
    }
    if (rho > errtol) {
        cerr << "\nGMRES failed to converge after " << maxIter_
             << " iterations.\n" << endl;
    }

    // y = H^-1 g.
    return arnoldi.qTimes(back_substitute(H, g, arnoldi.currentSize()));
}

/// \param[in,out] givensList Givens rotation series.
/// \param[in,out] H GMRES Hessenberg matrix.
/// \param[in] Hsize Integer, where \f$\mathbf{H}\f$ is (\c Hsize + 1) x \c
/// Hsize.
template <typename Vec>
void Gmres<Vec>::upperTriangulize(
    vector<GivensRotation>& givensList,
    MatrixXd& H,
    int Hsize
) const {
    int n{Hsize - 1};
    if (n > 0) {
        // Givens rotations applied to the last ("latest") column of H:
        for(int j = {0}; j < n ; ++j) {
            givensList[j].applyToACol(H, Hsize);
        }
    }

    double nu{std::sqrt(H(n, n) * H(n, n) + H(n + 1, n) * H(n + 1, n))};
    if (nu > 0) {
        givensList.push_back(
            GivensRotation{n, H(n, n) / nu, -H(n + 1, n) / nu}
        );
        assert(static_cast<int>(givensList.size()) == n + 1);
        givensList[n].applyToACol(H, Hsize);
    } else {
        // the new Givens rotation matrix is a 2 by 2 identity matrix.
        givensList.push_back(GivensRotation{n, 1, 0});
    }
}

/// \param[in] H Hessenberg matrix.
/// \param[in] g Equation right-hand side.
/// \param[in] Hsize Integer, where \f$\mathbf{H}\f$ is (\c Hsize + 1) x \c
/// Hsize.
/// \return An \c std::vector giving the solution \f$\mathbf{y}\f$.
template <typename Vec>
vector<double> Gmres<Vec>::back_substitute(
    const MatrixXd& H,
    const vector<double>& g,
    int Hsize
) const {
    int k{Hsize - 1}; // C++ index.
    assert(static_cast<int>(g.size()) == k + 2);
    vector<double> y(g.size() - 1, 0);

    for (int i = {k}; i >= 0; --i) {
        double sum{g[i]};
        for (int j = {i + 1}; j < k + 1; ++j) {
            sum -= H(i, j) * y[j];
        }
        y[i] = sum / H(i, i);
    }

    return y;
}

template <typename Vec>
inline void Gmres<Vec>::printConvergence(
    int iter,
    double rho,
    double errtol
) const {
    std::cout << "\nGMRES iteration " << iter << ": " << "||A * x - b|| = "
              << rho << "; tolerance = " << errtol << "." << std::endl;
}

/// Throw an \c InvalidValue exception if certain inputs don't make sense.
template <typename Vec>
void Gmres<Vec>::checkInputs() const {
    using aux_ns::InvalidValue;

    if (tol_ < 0) {
        throw InvalidValue{"tol >= 0 required."};
    }
    if (maxIter_ < 0) {
        throw InvalidValue{"maxIter >= 0 required."};
    }
}

// ************************************************************************** //
// GivensRotation class member functions.

/// The Givens rotation is the identity matrix, except for the 2 x 2 block
/// \f$\left[\begin{array}{cc} c & -s \\ s & c \end{array}\right]\f$ starting at
/// the \f$k\f$th row and \f$k\f$th column (with 1-indexing).
///
/// \param[in] k Index (with 1-indexing) where the rotation block is.
/// \param[in] c Upper left and lower right element of the 2 x 2 block.
/// \param[in] s Lower left element of the 2 x 2 block.
template <typename Vec>
Gmres<Vec>::GivensRotation::GivensRotation(int k, double c, double s) :
    k_{k},
    c_{c},
    s_{s}
{
}

/// Left-multiply the \f$i\f$th column (by 1-indexing) of the matrix
/// \f$\mathbf{H}\f$ by the Givens rotation matrix, and change the column
/// in-place.
///
/// \param [in,out] H Input matrix.
/// \param [in] i Column of \f$\mathbf{H}\f$ to apply the Givens rotation on.
/// \sa apply(vector<double>& v)
template <typename Vec>
inline void Gmres<Vec>::GivensRotation::applyToACol(MatrixXd& H, int i) {
    assert(H.rows() >= k_ + 2);
    twoDRotation(H(k_, i - 1), H(k_ + 1, i - 1));
}

/// This changes \f$\mathbf{v}\f$ in-place.
///
/// \param [in,out] v Input and output of the Givens rotation product.
/// \sa applyToACol(MatrixXd& H, int i)
template <typename Vec>
inline void Gmres<Vec>::GivensRotation::apply(vector<double>& v) {
    assert(static_cast<int>(v.size()) >= k_ + 2);
    twoDRotation(v[k_], v[k_+1]);
}
/**
/// The rotation applied is
/// \f{eqnarray*}{
///     v1 &\leftarrow& c * v1 - s * v2 \\
///     v2 &\leftarrow& s * v1 + c * v2
/// \f}
///
/// \param[in] v1 Input and output double.
/// \param[in] v2 Input and output double.
*/
template <typename Vec>
void Gmres<Vec>::GivensRotation::twoDRotation(
    double& v1,
    double& v2
) {
    double tmp = v1;

    v1 *= c_;
    v1 -= s_ * v2;

    v2 *= c_;
    v2 += s_ * tmp;
}

// ************************************************************************** //

} // End namespace linear_solver_ns

#endif // GMRES_H
