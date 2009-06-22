// EllipticSolver2d.cc
//
// Description:
// Implementation of the EllipticSolver2d class
//
// Author(s):
// Clancy Rowley
//
// Date: 17 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "EllipticSolver2d.h"
#include "VectorOperations.h"
#include <math.h>

namespace ibpm {
    
//------------------------------------------------------------------------------
// Elliptic solver (abstract base class)
//------------------------------------------------------------------------------

    EllipticSolver2d::EllipticSolver2d( int nx, int ny, double dx ) :
        // Need only interior points, so e-values are nx-1 by ny-1
        // 2D array of eigenvalues has indices that starts at 1
        _eigenvaluesOfInverse( nx-1, ny-1, 1, 1 ),
        _fft( nx-1, ny-1, 1, 1 ) {
        _nx = nx;
        _ny = ny;
        _dx = dx;
        _FFTWPlan = fftw_plan_r2r_2d( nx-1, ny-1, _fft, _fft,
            FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
    }
    
    EllipticSolver2d::~EllipticSolver2d() {
        fftw_destroy_plan( _FFTWPlan );
    }
    
    EllipticSolver2d::Array2d EllipticSolver2d::getLaplacianEigenvalues() const {
        // calculate eigenvalues of Laplacian
        double pi = 4. * atan(1.0);
        double bydx2 = 1. / (_dx * _dx);
        Array2d eig( _nx-1, _ny-1, 1, 1 );
        for (int i=1; i < _nx; ++i ) {
            for (int j=1; j < _ny; ++j ) {
                eig(i,j) = 2. * ( cos( (pi * i) / _nx ) + 
                                 cos( (pi * j) / _ny ) - 2. ) * bydx2;
            }
        }
        return eig;
    }
    
    // Take discrete sin transform of u, leaving result in v
    void EllipticSolver2d::sinTransform( const Array2d& u, Array2d& v ) const {
        assert( u.Size() == v.Size() );
        
        // copy input array to FFTW array
        for (unsigned int i = 0; i < u.Size(); ++i ) {
            _fft(i) = u(i);
        }
        
        fftw_execute( _FFTWPlan ); 
        
        // copy back to output array
        for (unsigned int i = 0; i < u.Size(); ++i) {
            v(i) = _fft(i);
        }
    }
    
    // Take inverse sin transform of u, leaving result in v
    // Note that inverse is the same as the forward transform, except for a
    // normalization factor
    void EllipticSolver2d::sinTransformInv( const Array2d& u, Array2d& v )
    const {
        assert( u.Size() == v.Size() );
        sinTransform( u, v );
        double normalizationFactor = 1. / (2 * _nx * 2 * _ny);
        for (unsigned int i=0; i < v.Size(); ++i ) {
            v(i) *= normalizationFactor;
        }
    }
    
    // Solve L u = f, single domain, assuming zero boundary conditions on u
    void EllipticSolver2d::solve(const Array2d& f, Array2d& u ) const {
        sinTransform( f, u );
        for (int i=1; i<_nx; ++i) {
            for (int j=1; j<_ny; ++j) {
                u(i,j) *= _eigenvaluesOfInverse(i,j);
            }
        }        
        sinTransformInv( u, u ); // normalize on inverse transform        
    }
        
    // Solve L u = f, with specified boundary conditions on u.
    // Note that u contains only the interior points of the domain.
    void EllipticSolver2d::solve( const Array2d& f, const BC& bc, Array2d& u )
        const {
        Array2d& rhs = u;  // use u as storage for rhs of Poisson equation
        getRHS( f, bc, rhs );
        solve( rhs, u );
    }
    
//------------------------------------------------------------------------------
// Poisson solver
//------------------------------------------------------------------------------    
    
    PoissonSolver2d::PoissonSolver2d( int nx, int ny, double dx ) :
        EllipticSolver2d( nx, ny, dx ) {
        calculateEigenvalues();
    }
    
    void PoissonSolver2d::calculateEigenvalues() {
        Array2d eigL = getLaplacianEigenvalues();
        for (int i=1; i < _nx; ++i ) {
            for (int j=1; j < _ny; ++j ) {
                _eigenvaluesOfInverse(i,j) = 1. / eigL(i,j);
            }
        }
    }
        
    // Set rhs = f - L * bc, in preparation for Poisson solve
    void PoissonSolver2d::getRHS( const Array2d& f, const BC& bc, Array2d& rhs ) const {
        assert( f.Nx() == rhs.Nx() );
        assert( f.Ny() == rhs.Ny() );
        // if input and output arrays are not the same, copy data to rhs
        if ( &rhs != &f ) {
            for (unsigned int i=0; i<f.Size(); ++i) {
                rhs(i) = f(i);
            }
        }
        // subtract L(bc) from rhs
        const double byDx2 = 1 / (_dx * _dx);
        for (int i=1; i<_nx; ++i) {
            rhs(i,1) -= bc.bottom(i) * byDx2;
            rhs(i,_ny-1) -= bc.top(i) * byDx2;
        }
        for (int j=1; j<_ny; ++j) {
            rhs(1,j) -= bc.left(j) * byDx2;
            rhs(_nx-1,j) -= bc.right(j) * byDx2;
        }
    }
    
//------------------------------------------------------------------------------
// Helmholtz solver
//------------------------------------------------------------------------------    
    
    HelmholtzSolver2d::HelmholtzSolver2d( int nx, int ny, double dx, double alpha ) :
        EllipticSolver2d( nx, ny, dx ),
        _alpha(alpha) {
        calculateEigenvalues();
    }
    
    void HelmholtzSolver2d::calculateEigenvalues() {
        Array2d eigL = getLaplacianEigenvalues();
        for (int i=1; i < _nx; ++i ) {
            for (int j=1; j < _ny; ++j ) {
                _eigenvaluesOfInverse(i,j) = 1. / (1 + _alpha * eigL(i,j));
            }
        }        
    }

    // Set rhs = f - alpha * L * bc, in preparation for Helmholtz solve
    void HelmholtzSolver2d::getRHS( const Array2d& f, const BC& bc, Array2d& rhs ) const {
        assert( f.Nx() == rhs.Nx() );
        assert( f.Ny() == rhs.Ny() );
        // if input and output arrays are not the same, copy data to rhs
        if ( &rhs != &f ) {
            for (unsigned int i=0; i<f.Size(); ++i) {
                rhs(i) = f(i);
            }
        }
        // subtract alpha * L(bc) from rhs
        const double alphaByDx2 = _alpha / (_dx * _dx);
        for (int i=1; i<_nx; ++i) {
            rhs(i,1) -= bc.bottom(i) * alphaByDx2;
            rhs(i,_ny-1) -= bc.top(i) * alphaByDx2;
        }
        for (int j=1; j<_ny; ++j) {
            rhs(1,j) -= bc.left(j) * alphaByDx2;
            rhs(_nx-1,j) -= bc.right(j) * alphaByDx2;
        }
    }
    
    
} // namespace ibpm
