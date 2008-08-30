// CholeskySolver.cc
//
// Description:
// Implementation of the CholeskySolver class
//
// Author(s):
// Clancy Rowley
//
// Date: 28 Aug 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "CholeskySolver.h"
#include "NavierStokesModel.h"
#include "assert.h"
#include <iostream>

using namespace std;

namespace ibpm {

// Allocate memory for the Cholesky factorization, but do not compute it
CholeskySolver::CholeskySolver(
    const NavierStokesModel& model,
    double alpha
    ) :
    ProjectionSolver( model, alpha ),
    _numPoints( model.getGeometry()->getNumPoints() ),
    _size( 2 * _numPoints ),
    _lower( _size ),
    _diagonal( _size ) {
    computeCholesky();
}

CholeskySolver::~CholeskySolver() {
    // Deallocate memory for the Cholesky factorization
    // (no need when Blitz++ arrays are initialized as above)
}

// Compute the matrix M to be factored
void CholeskySolver::computeMatrixM( Array<double,2>& matrixM ) {
    BoundaryVector e( _numPoints ); // basis vector
    BoundaryVector x( _numPoints ); // x = M(e)
        
    for ( int j=0; j<_size; ++j ) {
        // Compute j-th column of M
        e = 0;
        e(j) = 1;       // j-th basis vector
        M( e, x );      // Compute x = M(e)
        // Copy into matrix M
        for ( int i=0; i<_size; ++i ) {
            matrixM(i,j) = x(i);
        }
    }
}

// 
void CholeskySolver::computeFactorization( const Array<double,2>& matrixM ) {
    
    _lower = matrixM;
    for ( int i=0; i<_size; ++i ) {
        for ( int j=0; j<_size; ++j ) {
            double sum = _lower(i,j);
            for ( int k=i-1; k>=0; --k ) {
                sum -= _lower(i,k) * _lower(j,k);
            }
            if( i==j ) {
                assert( sum > 0 );
                _diagonal(i) = sqrt(sum);
            }
            else {
                _lower(j,i) = sum / _diagonal(i);
            }
        }
    }
}

void CholeskySolver::computeCholesky() {
    // Build matrix M explicitly, one column at a time
    Array<double,2> matrixM( _size );
    computeMatrixM( matrixM );

    // Compute Cholesky factorization
    computeFactorization( matrixM );
}

// Load a Cholesky decomposition from the specified file.
// Return true if successful
bool CholeskySolver::loadCholesky(const string& filename) {
    // Not implemented
    return false;
}

// Save a Cholesky decomposition to the specified file, overwriting if
// necessary.
// Return true if successful
bool CholeskySolver::saveCholesky(const string& filename) {
    // Not implemented
    return false;
}

// Solve A x = b using the Cholesky factorization A = LL*
void CholeskySolver::Minv(
    const BoundaryVector& b,
    BoundaryVector& x
    ) {

    // Solve L y = b for y
    // (Here, y and x use the same memory, for efficiency)
    for ( int i=0; i<_size; ++i ) {
        double sum = b(i);
        for ( int k=i-1; k>=0; --k ) {
            sum -= _lower(i,k) * x(k);
        }
        x(i) = sum / _diagonal(i);
    }

    // Solve L^Tx = y for x
    for ( int i=_size-1; i>=0; --i ) {
        double sum = x(i);
        for (int k=i+1; k<_size; ++k ) {
            sum -= _lower(k,i) * x(k);
        }
        x(i) = sum / _diagonal(i);
    }
}

} // namespace ibpm