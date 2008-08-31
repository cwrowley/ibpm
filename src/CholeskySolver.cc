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
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

namespace ibpm {

// Allocate memory for the Cholesky factorization, but do not compute it
CholeskySolver::CholeskySolver(
    const NavierStokesModel& model,
    double alpha
    ) :
    ProjectionSolver( model, alpha ),
    _numPoints( model.getGeometry().getNumPoints() ),
    _size( 2 * _numPoints ),
    _lower( _size ),
    _diagonal( _size ),
    _hasBeenInitialized( false )
{}

CholeskySolver::~CholeskySolver() {
    // Deallocate memory for the Cholesky factorization
    // (no need when Blitz++ arrays are initialized as above)
}

void CholeskySolver::init() {
    // Build matrix M explicitly, one column at a time
    Array<double,2> matrixM( _size );
    computeMatrixM( matrixM );

    // Compute Cholesky factorization
    computeFactorization( matrixM );
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

// Compute the Cholesky factorization of matrixM
//    M = L L*
// where L is lower triangular.
// Preconditions:
//      M is symmetric
// Postconditions:
//      _lower contains the strictly lower triangular part of L (no diagonal)
//      _diag  contains the diagonal elements of L
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
    _hasBeenInitialized = true;
}

// Load a Cholesky decomposition from the specified file.
// Return true if successful
bool CholeskySolver::loadCholesky(const string& filename) {
    ifstream infile( filename.c_str() );
    if ( ! infile.good() ) return false;
    int n;
    infile >> n;
    if (n != _size) return false;
    for ( int i=0; i<_size; ++i ) {
        infile >> _diagonal(i);
    }

    // check the marker, to make sure we did not get off track
    char c;
    infile >> c;
    if (c != '#') return false;

    // read the lower triangular portion
    for ( int i=0; i<_size; ++i) {
        for (int j=0; j<i; ++j) {
            infile >> _lower(i,j);
        }
    }
    _hasBeenInitialized = true;
    return true;
}

// Save a Cholesky decomposition to the specified file, overwriting if
// necessary.
// Return true if successful
// Note: saves only strictly lower triangular portion of _lower, since
// that is all that is needed for back substitution
bool CholeskySolver::saveCholesky(const string& filename) {
    assert( _hasBeenInitialized );
    ofstream outfile( filename.c_str() );
    if ( ! outfile.good() ) return false;
    outfile << _size << endl;
    // write the diagonal part
    for ( int i=0; i<_size; ++i ) {
        outfile << setprecision(17) << _diagonal(i) << endl;
    }

    // insert marker, as a crude verification that we did not skip a line
    // when reading back in
    outfile << "#" << endl;

    // write the lower triangular part
    for ( int i=0; i<_size; ++i) {
        for (int j=0; j<i; ++j) {
            outfile << setprecision(17) << _lower(i,j) << endl;
        }
    }
    return true;
}

// Solve A x = b using the Cholesky factorization A = LL*
void CholeskySolver::Minv(
    const BoundaryVector& b,
    BoundaryVector& x
    ) {

    assert( _hasBeenInitialized );
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