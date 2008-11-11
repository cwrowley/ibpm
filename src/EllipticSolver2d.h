#ifndef _ELLIPTICSOLVER_H_
#define _ELLIPTICSOLVER_H_

#include "Array.h"
#include "BC.h"
#include <fftw3.h>

namespace ibpm {
    
/*!
    \file EllipticSolver2d.h
    \class EllipticSolver2d
     
    \brief Class for solving Poisson and Helmholtz equations on a uniform grid
     
    Uses a sin transform, or a nested family of sin transforms.  Note that this
    is an abstract base class, and cannot be instantiated.

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 17 Sep 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class EllipticSolver2d {
public:
    /// Type for arrays used to store 2d scalar fields
    typedef Array::Array2<double> Array2d;
    
    /// \brief Instantiate a new solver.
    ///
    /// Here, nx and ny are the number of cells in x- and y-directions, and dx
    /// is the grid spacing.
    /// Note that the number of gridpoints is (nx+1) x (ny+1), and the number
    /// of interior gridpoints is (nx-1) x (ny-1).
    EllipticSolver2d( int nx, int ny, double dx );
    virtual ~EllipticSolver2d() = 0;

    /// \brief Solve L u = f, assuming zero boundary conditions on u
    /// The 2D arrays f and u must have indices (1..nx-1, 1..ny-1)
    void solve( const Array2d& f, Array2d& u ) const;

    /// \brief Solve L u = f, with specified boundary conditions on u.
    /// Note that u contains only the interior points of the domain.
    /// The 2D arrays f and u must have indices (1..nx-1, 1..ny-1)
    /// and the BC object has size (nx,ny)
    void solve( const Array2d& f, const BC& bc, Array2d& u ) const;

protected:
    Array2d getLaplacianEigenvalues() const;
    virtual void getRHS( const Array2d& f, const BC& bc, Array2d& rhs ) const = 0;
    Array2d _eigenvaluesOfInverse;
    int _nx;
    int _ny;
    double _dx;

private:
    void sinTransform( const Array2d& u, Array2d& v ) const;
    void sinTransformInv( const Array2d& u, Array2d& v ) const;
    fftw_plan _FFTWPlan;
    Array2d _fft;
};

/******************************************************************************/

/*! \class PoissonSolver2d
    
    \brief Solve a Poisson equation on a uniform grid.
    Solves L u = f, with specified boundary conditions on u, where
    L is the Laplacian.
*/
class PoissonSolver2d : public EllipticSolver2d {
public:
    PoissonSolver2d( int nx, int ny, double dx );
protected:
    // set rhs = f - L * bc
    void getRHS( const Array2d& f, const BC& bc, Array2d& rhs ) const;
private:
    void calculateEigenvalues();
};

/******************************************************************************/
    
/*! \class HelmholtzSolver2d
    \brief Solve a Helmholtz equation on a uniform grid
    Solves (1 + alpha * L) u = f, with specified boundary conditions on u,
    where L is the Laplacian.
*/
class HelmholtzSolver2d : public EllipticSolver2d {
public:
    HelmholtzSolver2d( int nx, int ny, double dx, double alpha );
protected:
    // set rhs = f - alpha * L * bc
    void getRHS( const Array2d& f, const BC& bc, Array2d& rhs ) const;
private:
    double _alpha;
    void calculateEigenvalues();
};

} // namespace ibpm

#endif /* _ELLIPTICSOLVER_H_ */
