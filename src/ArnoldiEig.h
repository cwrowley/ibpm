/// \file ArnoldiEig.h
/// \brief A template for estimating an eigendecomposition by Arnoldi
/// decomposition.
/// \author Kevin K. Chen, Jonathan H. Tu (Princeton University)

#ifndef ARNOLDI_EIG_H
#define ARNOLDI_EIG_H

#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <Eigen>
#include <Eigenvalues>

#include "ArnoldiDecomp.h"
#include "exceptions.h"
#include "templateAliases.h"

// ************************************************************************** //

/// \namespace eig_ns
/// \brief The namespace used for all eigenvalue-related constructs.
namespace eig_ns {

using std::unique_ptr;
using std::vector;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::VectorXcd;

using arnoldi_ns::ArnoldiDecomp;

typedef std::complex<double> cdouble; ///< A complex double.

// ************************************************************************** //

/// \brief Indicator for whether the eigenvalues should be sorted by real part
/// or absolute value.
enum class SortByType {real, abs};

// ************************************************************************** //
// ArnoldiEig template declaration.

/// A template for estimating an eigendecomposition by Arnoldi decomposition.

/// This template interacts with an ArnoldiDecomp object to find the approximate
/// eigenvalues and eigenvectors of an operator \f$\mathbf{A}\f$.  Any number of
/// eigenvalues and eigenvectors, up to the size of the Arnoldi Hessenberg
/// matrix, can be computed.
///
/// \tparam Vec Vector type for \c ArnoldiEig.
/// \sa Eig, EmpiricalEig
template <typename Vec>
class ArnoldiEig {
public:
    /// Constructor.
    ArnoldiEig(
        const alias::Function<Vec>& f,
        const alias::InnerProduct<Vec>& ip,
        const Vec& b,
        int n,
        SortByType sort_by = SortByType::real,
        bool verbose = false
    );
    ~ArnoldiEig() noexcept {} ///< Destructor (empty).

    VectorXcd eVals(); ///< \brief Compute eigenvalues.
    /// Compute eigenvalues and eigenvectors.
    void eig(VectorXcd& evalues, vector<vector<Vec>>& evectors);
    /// Compute eigenvalues and and a subset of eigenvectors.
    void eig(VectorXcd& evalues, vector<vector<Vec>>& evectors, int n_vectors);

    /// Sort eigenvalues and eigenvectors by eigenvalue real part or magnitude.
    void sortEigPair(VectorXcd& evalues, MatrixXcd& evectors) const;

private:
    class EigPair; ///< A class holding one eigenvalue and its eigenvector.

    ArnoldiDecomp<Vec> arnoldi_decomp_; ///< Arnoldi decomposition object.
    int n_; ///< Number of eigenvalues to compute.
    const Vec v_; ///< Any \c Vec object, for copying.
    const SortByType sort_by_; ///< \c real or \c abs; how to sort eigenvalues.

    ArnoldiEig(const ArnoldiEig&) = delete;
    ArnoldiEig(ArnoldiEig&&) = delete;

    void operator=(const ArnoldiEig&) = delete;
    void operator=(ArnoldiEig&&) = delete;

    void checkInputs() const; ///< Check for valid inputs.

    /// Get the upper left corner of \f$\mathbf{H}_m\f$.
    MatrixXd getH(int m);

    /// \brief Make evectors have the right size for
    /// eig(VectorXcd& evalues, vector<vector<Vec>>& evectors).
    void initEVectors(vector<vector<Vec>>& evectors, int n) const;
    /// Sort a complex vector by descending real part or magnitude.
    void sortVectorXcd(VectorXcd& x) const;

    /// Compare complex numbers' real parts.
    static bool realLessThan(const cdouble& x, const cdouble& y);
    /// Compare complex numbers' magnitudes.
    static bool absLessThan(const cdouble& x, const cdouble& y);
};

// ************************************************************************** //
// EigPair class declaration.

/// The \c EigPair class is used to sort eigenvalue and eigenvector pairs by the
/// value of the eigenvalue, which can be done by real part or absolute value.
///
/// \tparam Vec Vector type for \c EigPair.
template <typename Vec>
class ArnoldiEig<Vec>::EigPair {
public:
    /// Make an eigenvalue-eigenvector pair.
    EigPair(const cdouble& eValue, const VectorXcd& eVector);
    EigPair(EigPair&&) = default; ///< Move constructor.

    ~EigPair() noexcept {} ///< Destructor (empty).

    /// Set an <tt>EigPair</tt>'s eigenvalue and eigenvector to another's.
    EigPair& operator=(EigPair&&) = default;

    const cdouble& getEValue() const; ///< Get the eigenvalue from the pair.
    const VectorXcd& getEVector() const; ///< Get the eigenvector from the pair.

    /// Compare the real parts of the eigenvalues of two <tt>EigPair</tt>s.
    static bool realLessThan(
        const unique_ptr<EigPair>& ep1,
        const unique_ptr<EigPair>& ep2
    );
    /// Compare the magnitudes of the eigenvalues of two <tt>EigPair</tt>s.
    static bool absLessThan(
        const unique_ptr<EigPair>& ep1,
        const unique_ptr<EigPair>& ep2
    );

private:
    const cdouble eValue_; ///< The eigenvalue.
    const VectorXcd eVector_; ///< The eigenvector.

    EigPair(const EigPair&) = delete;

    EigPair& operator=(const EigPair&) = delete;
};

// ************************************************************************** //
// ArnoldiEig template member functions.

/// \param[in] f The \c std::function to perform an eigendecomposition on.
/// \param[in] ip The \c std::function defining the inner product for \c Vec.
/// \param[in] b The initial guess for the Arnoldi decomposition.
/// \param[in] n The number of eigenvalues to estimate.
/// \param[in] sort_by \c SortByType::real or \c SortByType::abs, defining how
/// to sort the eigenvalues and their corresponding eigenvectors.  Use \c
/// SortByType::real for continuous-time dynamics and \c SortByType::abs for
/// discrete-time dynamics.
/// \param[in] verbose \c true to activate verbose printing.
template <typename Vec>
ArnoldiEig<Vec>::ArnoldiEig(
    const alias::Function<Vec>& f,
    const alias::InnerProduct<Vec>& ip,
    const Vec& b,
    int n,
    SortByType sort_by,
    bool verbose
) :
    arnoldi_decomp_{f, ip, b, n, verbose},
    n_{n},
    v_{b},
    sort_by_{sort_by}
{
}

/// This does not compute eigenvectors, and is faster than <tt>eig(VectorXcd&
/// evalues, vector<vector<Vec>>& evectors)</tt>.
///
/// \return Complex vector of eigenvalues.
/// \sa eig(VectorXcd& evalues, vector<vector<Vec>>& evectors),
/// eig(VectorXcd& evalues, vector<vector<Vec>>& evectors, int n)
template <typename Vec>
VectorXcd ArnoldiEig<Vec>::eVals() {
    MatrixXd H{getH(n_)};
    Eigen::EigenSolver<MatrixXd> eigenSolver{H, false};

    VectorXcd eigenvalues{eigenSolver.eigenvalues()};
    sortVectorXcd(eigenvalues); // Sort by descending real part.

    return eigenvalues;
}

/// \param[out] evalues Complex vector of eigenvalues.
/// \param[out] evectors Vector of vector of eigenvectors.
/// <tt>evectors[i][0]</tt> is the real part of the \f$i\f$th eigenvector;
/// <tt>evectors[i][1]</tt> is the imaginary part.
///
/// \sa eVals(), eig(VectorXcd& evalues, vector<vector<Vec>>& evectors, int n)
template <typename Vec>
inline void ArnoldiEig<Vec>::eig(
    VectorXcd& evalues,
    vector<vector<Vec>>& evectors
) {
    eig(evalues, evectors, n_);
}

/// \param[out] evalues Complex vector of eigenvalues.
/// \param[out] evectors Vector of vector of eigenvectors.
/// <tt>evectors[i][0]</tt> is the real part of the \f$i\f$th eigenvector;
/// <tt>evectors[i][1]</tt> is the imaginary part.
/// \param[in] n_vectors Number of eigenvectors to return.  Regardless of \c
/// n_vectors, the number of eigenvalues to solve will always be the value of \c
/// n given to the class constructor.
///
/// \sa eVals(), eig(VectorXcd& evalues, vector<vector<Vec>>& evectors)
template <typename Vec>
void ArnoldiEig<Vec>::eig(
    VectorXcd& evalues,
    vector<vector<Vec>>& evectors,
    int n_vectors
) {
    if (n_vectors < 0 || n_vectors > n_) {
        std::ostringstream s{};
        s << "n_vectors = " << n_vectors << " must be in the range [0, " << n_
          << "].";
        throw aux_ns::InvalidValue(s.str());
    }

    // Initialize evectors with the right size.
    initEVectors(evectors, n_vectors);

    // Eigen::EigenSolver<MatrixXd> cannot be constructed with eigenvector
    // solving when the matrix is empty.  Just quit early instead.
    if (n_ == 0) {
        VectorXcd empty{};
        evalues = empty;
        return;
    }

    MatrixXd H{getH(n_)};
    Eigen::EigenSolver<MatrixXd> eigenSolver{H, true};

    evalues = eigenSolver.eigenvalues();
    MatrixXcd V{eigenSolver.eigenvectors()};

    // Sort eigenvalues and eigenvectors by descending eigenvalue real part.
    sortEigPair(evalues, V);

    for (int i = {0}; i < n_vectors; ++i) {
        evectors[i] = arnoldi_decomp_.qTimes(static_cast<VectorXcd>(V.col(i)));
    }
}

/// \param[in] evalues Vector of eigenvalues.
/// \param[in] evectors Matrix of eigenvectors.
template <typename Vec>
void ArnoldiEig<Vec>::sortEigPair(
    VectorXcd& evalues,
    MatrixXcd& evectors
) const {
    using std::sort;

    int n{static_cast<int>(evalues.size())};

    // Check that sizes make sense.
    if (evectors.rows() != evectors.cols()) {
        std::cout << "Warning: evectors is not square." << std::endl;
    }
    if (n != evectors.cols()) {
        std::ostringstream s{};
        s << "Error: evalues must have as many elements (" << n
          << ") as evectors has columns (" << evectors.cols() << ").";
        throw std::length_error{s.str()};
    }

    // Initialize pairs with n copies null EigPair pointers.
    vector<unique_ptr<EigPair>> pairs(n);
    // Fill up pairs.
    for (int j = {0}; j < n; ++j) {
        pairs[j].reset(
            new EigPair{
                evalues.coeff(j),
                static_cast<VectorXcd>(evectors.col(j))
            }
        );
    }

    // Sort pairs by eigenvalue real part or magnitude.
    if (sort_by_ == SortByType::real) {
        sort(pairs.begin(), pairs.end(), EigPair::realLessThan);
    } else { // sort_by_ is SortByType::abs.
        sort(pairs.begin(), pairs.end(), EigPair::absLessThan);
    }
    std::reverse(pairs.begin(), pairs.end());

    // Set evalues and evectors to the sorted values in pairs.
    for (int j = {0}; j < n; ++j) {
        evalues[j] = pairs[j]->getEValue();
        evectors.col(j) = pairs[j]->getEVector();
    }
}

/// Throw an \c InvalidValue exception if certain inputs don't make sense.
template <typename Vec>
void ArnoldiEig<Vec>::checkInputs() const {
    if (n_ < 0) {
        std::ostringstream s{};
        s << "n = " << n_ << " must be non-negative.";
        throw aux_ns::InvalidValue{s.str()};
    }
}

/// Iterate the Arnoldi decomposition until it has a size of \f$m\f$, and then
/// return the upper left \f$m \times m\f$ corner of the Hessenberg matrix.
///
/// \f$\mathbf{H}_m\f$.
template <typename Vec>
inline MatrixXd ArnoldiEig<Vec>::getH(int m) {
    arnoldi_decomp_.stepTo(m);
    return arnoldi_decomp_.getH().topLeftCorner(m, m);
}

/// \param[in] evectors The eigenvectors object.
/// \param[in] n The number of eigenvectors to solve.
template <typename Vec>
void ArnoldiEig<Vec>::initEVectors(vector<vector<Vec>>& evectors, int n) const
{
    evectors.erase(evectors.begin(), evectors.end());
    evectors.insert(evectors.end(), n, vector<Vec>(2, v_));
}

/// The vector is sorted either by reverse real part of reverse magnitude,
/// depending on the value of \c sort_by_.
///
/// \param[in,out] x The vector to sort.
template <typename Vec>
void ArnoldiEig<Vec>::sortVectorXcd(VectorXcd& x) const {
    using std::sort;

    if (sort_by_ == SortByType::real) {
        sort(x.data(), x.data() + x.size(), ArnoldiEig<Vec>::realLessThan);
    } else { // sort_by_ is SortByType::abs.
        sort(x.data(), x.data() + x.size(), ArnoldiEig<Vec>::absLessThan);
    }

    std::reverse(x.data(), x.data() + x.size());
}

/// \param[in] x A complex scalar.
/// \param[in] y A complex scalar.
///
/// \return True if \f$\mathrm{Re} (x) < \mathrm{Re} (y)\f$; false otherwise.
/// \sa absLessThan
template <typename Vec>
inline bool ArnoldiEig<Vec>::realLessThan(const cdouble& x, const cdouble& y) {
    return x.real() < y.real();
}

/// \param[in] x A complex scalar.
/// \param[in] y A complex scalar.
///
/// \return True if \f$|x| < |y|\f$; false otherwise.
/// \sa realLessThan
template <typename Vec>
inline bool ArnoldiEig<Vec>::absLessThan(const cdouble& x, const cdouble& y) {
    using std::abs;

    return abs(x) < abs(y);
}


// ************************************************************************** //
// EigPair class member functions.

/// \param[in] eValue the eigenvalue.
/// \param[in] eVector the eigenvector.
template <typename Vec>
inline ArnoldiEig<Vec>::EigPair::EigPair(
    const cdouble& eValue,
    const VectorXcd& eVector
) :
    eValue_(eValue),
    eVector_(eVector)
{
}

/// \return The eigenvalue.
/// \sa getEValue()
template <typename Vec>
inline const cdouble& ArnoldiEig<Vec>::EigPair::getEValue() const {
    return eValue_;
}

/// \return The eigenvector.
/// \sa getEVector()
template <typename Vec>
inline const VectorXcd& ArnoldiEig<Vec>::EigPair::getEVector() const {
    return eVector_;
}

/// \return \c true if <tt>ep1</tt>'s eigenvalue real part is less than
/// <tt>ep2</tt>'s; \c false otherwise.
/// \sa absLessThan(const EigPair& ep1, const EigPair& ep2)
template <typename Vec>
inline bool ArnoldiEig<Vec>::EigPair::realLessThan(
    const unique_ptr<EigPair>& ep1,
    const unique_ptr<EigPair>& ep2
) {
    return ArnoldiEig<Vec>::realLessThan(ep1->getEValue(), ep2->getEValue());
}

/// \return \c true if <tt>ep1</tt>'s eigenvalue magnitude is less than
/// <tt>ep2</tt>'s; \c false otherwise.
/// \sa realLessThan(const EigPair& ep1, const EigPair& ep2)
template <typename Vec>
inline bool ArnoldiEig<Vec>::EigPair::absLessThan(
    const unique_ptr<EigPair>& ep1,
    const unique_ptr<EigPair>& ep2
) {
    return ArnoldiEig<Vec>::absLessThan(ep1->getEValue(), ep2->getEValue());
}

// ************************************************************************** //

} // End namespace arnoldi_ns

#endif // ARNOLDI_EIG_H
