// EllipticSolver.cc
//
// Description:
// Implementation of the EllipticSolver class
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

#include "EllipticSolver.h"
#include "VectorOperations.h"

using Array::Array2;

namespace ibpm {
    EllipticSolver::EllipticSolver( int nx, int ny, double dx ) :
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
    
    EllipticSolver::EllipticSolver( const Grid& grid ) :
        // Need only interior points, so e-values are nx-1 by ny-1
        // 2D array of eigenvalues has indices that starts at 1
        _eigenvaluesOfInverse( grid.Nx()-1, grid.Ny()-1, 1, 1 ),
        _fft( grid.Nx()-1, grid.Ny()-1, 1, 1 ) {
        _nx = grid.Nx();
        _ny = grid.Ny();
        _dx = grid.Dx();
        _FFTWPlan = fftw_plan_r2r_2d( grid.Nx()-1, grid.Ny()-1, _fft, _fft,
            FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
    }
    
    EllipticSolver::~EllipticSolver() {
        fftw_destroy_plan( _FFTWPlan );
    }
    
    // Solve L u = f, assuming zero boundary conditions on u
    void EllipticSolver::solve( const Scalar& f, Scalar& u ) const {
        Scalar v( u.getGrid() );
        sinTransform( f, v );
        for (int i=1; i<_nx; ++i) {
            for (int j=1; j<_ny; ++j) {
                 v(i,j) *= _eigenvaluesOfInverse(i,j);
            }
        }        
        sinTransformInv( v, u ); // normalize on inverse transform
    }
    
    Scalar EllipticSolver::solve( const Scalar& f ) const {
        Scalar u( f.getGrid() );
        solve( f, u );
        return u;
    }

    Array2<double> EllipticSolver::getLaplacianEigenvalues() const {
        // calculate eigenvalues of Laplacian
        double pi = 4. * atan(1.0);
        double bydx2 = 1. / (_dx * _dx);
        Array2<double> eig( _nx-1, _ny-1, 1, 1 );
        for (int i=1; i < _nx; ++i ) {
            for (int j=1; j < _ny; ++j ) {
                eig(i,j) = 2. * ( cos( (pi * i) / _nx ) + 
                                  cos( (pi * j) / _ny ) - 2. ) * bydx2;
            }
        }
        return eig;
    }
    
    // Solve L u = f, with specified boundary conditions on u.
    // Note that u contains only the interior points of the domain.
//    void EllipticSolver::solve( const Scalar& f, const BC& bc, Scalar& u ) {
//        Scalar rhs( f.getGrid() );
//        getRHS( f, bc, rhs );
//        solve( rhs, u );
//    }
    
    void EllipticSolver::sinTransform( const Scalar& u, Scalar& v ) const {        
        // copy input array to FFTW array
        for (int i = 1; i < _nx; ++i) {
            for (int j = 1; j < _ny; ++j) {
                _fft(i,j) = u(i,j);
            }
        }
        
        fftw_execute( _FFTWPlan ); 
        
        // copy back to output array
        for (int i = 1; i < _nx; ++i) {
            for (int j = 1; j < _ny; ++j) {
                v(i,j) = _fft(i,j);
            }     
        }
        // set boundary terms to zero
        for (int i=0; i<= _nx; ++i) {
            v(i,0) = 0.;
            v(i,_ny) = 0.;
        }
        for (int j=0; j <= _ny; ++j) {
            v(0,j) = 0.;
            v(_nx,j) = 0.;
        }
    }

    void EllipticSolver::sinTransformInv( const Scalar& u, Scalar& v ) const {
        sinTransform( u, v );
        double normalizationFactor = 1. / (2 * _nx * 2 * _ny);
        v *= normalizationFactor;
    }
    
    void PoissonSolver::calculateEigenvalues() {
        Array2<double> eigL = getLaplacianEigenvalues();
        for (int i=1; i < _nx; ++i ) {
            for (int j=1; j < _ny; ++j ) {
                _eigenvaluesOfInverse(i,j) = 1. / eigL(i,j);
            }
        }
    }
        
    PoissonSolver::PoissonSolver( int nx, int ny, double dx ) :
        EllipticSolver( nx, ny, dx ) {
        calculateEigenvalues();
    }

    PoissonSolver::PoissonSolver( const Grid& grid ) :
        EllipticSolver( grid ) {
        calculateEigenvalues();
    }
    
//    void PoissonSolver::getRHS( const Scalar& f, const BC& bc, Scalar& rhs ) {
//        // TODO: set rhs to f - L * bc
//    }
    
    void HelmholtzSolver::calculateEigenvalues() {
        Array2<double> eigL = getLaplacianEigenvalues();
        for (int i=1; i < _nx; ++i ) {
            for (int j=1; j < _ny; ++j ) {
                _eigenvaluesOfInverse(i,j) = 1. / (1 + _alpha * eigL(i,j));
            }
        }        
    }

    HelmholtzSolver::HelmholtzSolver( const Grid& grid, double alpha ) :
        EllipticSolver( grid ),
        _alpha(alpha) {
        calculateEigenvalues();
    }

    HelmholtzSolver::HelmholtzSolver( int nx, int ny, double dx, double alpha ) :
        EllipticSolver( nx, ny, dx ),
        _alpha(alpha) {
        calculateEigenvalues();
    }

    
//    void HelmholtzSolver::getRHS( const Scalar& f, const BC& bc, Scalar& rhs ) {
//        // TODO: set rhs to f - alpha L * bc
//    }

    
} // namespace ibpm