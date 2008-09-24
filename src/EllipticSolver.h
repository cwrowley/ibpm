#ifndef _ELLIPTICSOLVER_H_
#define _ELLIPTICSOLVER_H_

#include "Array.h"
#include "Grid.h"
#include "Scalar.h"
#include <fftw3.h>
//#include "BC.h"

using Array::Array2;

namespace ibpm {
    
/*!
    \file EllipticSolver.h
    \class EllipticSolver
     
    \brief Class for solving Poisson and Helmholtz equations on a uniform grid
     
    Uses a sin transform, or a nested family of sin transforms.  Note that this
    is an abstract base class, and cannot be instantiated.

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 17 Sep 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class EllipticSolver {
public:
    /// \brief Instantiate a new solver.
    ///
    /// Here, nx and ny are the number of cells in x- and y-directions.  Note
    /// that the number of gridpoints is (nx+1) x (ny+1), and the number of
    /// interior gridpoints is (nx-1) x (ny-1).
    EllipticSolver( int nx, int ny, double dx );
    EllipticSolver( const Grid& grid );
    virtual ~EllipticSolver() = 0;

    /// \brief Solve L u = f, assuming zero boundary conditions on u
    void solve( const Scalar& f, Scalar& u ) const;
    Scalar solve( const Scalar& f ) const;

    /// \brief Solve L u = f, with specified boundary conditions on u.
    /// Note that u contains only the interior points of the domain.
//    void solve( const Scalar& f, const BC& bc, Scalar& u );

protected:
    Array2<double> getLaplacianEigenvalues() const;
    //    void getRHS( const Scalar& f, const BC& bc, Scalar& rhs ) = 0;
    Array2<double> _eigenvaluesOfInverse;
    int _nx;
    int _ny;
    double _dx;

private:
    void EllipticSolver::sinTransform( const Scalar& u, Scalar& v ) const;
    void EllipticSolver::sinTransformInv( const Scalar& u, Scalar& v ) const;
    fftw_plan _FFTWPlan;
    Array2<double> _fft;
};

/******************************************************************************/

/*! \class PoissonSolver
    
    \brief Solve a Poisson equation on a uniform grid.
    Solves L u = f, with specified boundary conditions on u, where
    L is the Laplacian.
*/
class PoissonSolver : public EllipticSolver {
public:
    PoissonSolver( int nx, int ny, double dx );
    PoissonSolver( const Grid& grid );
//    getRHS( const Scalar& f, const BC& bc, Scalar& rhs );
private:
    void calculateEigenvalues();
};

/******************************************************************************/
    
/*! \class HelmholtzSolver
    \brief Solve a Helmholtz equation on a uniform grid
    Solves (1 + alpha * L) u = f, with specified boundary conditions on u,
    where L is the Laplacian.
*/
class HelmholtzSolver : public EllipticSolver {
public:
    HelmholtzSolver( int nx, int ny, double dx, double alpha );
    HelmholtzSolver( const Grid& grid, double alpha );
//    getRHS( const Scalar& f, const BC& bc, Scalar& rhs );
private:
    double _alpha;
    void calculateEigenvalues();
};

} // namespace ibpm

#endif /* _ELLIPTICSOLVER_H_ */
