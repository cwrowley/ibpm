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

using std::ifstream;
using std::ofstream;
using std::setprecision;

namespace ibpm {

// Allocate memory for the Cholesky factorization, but do not compute it
CholeskySolver::CholeskySolver(
    const Grid& grid,
    const NavierStokesModel& model,
    double beta
    ) :
    ProjectionSolver( grid, model, beta ),
    _numPoints( model.getNumPoints() ),
    _size( 2 * _numPoints ),
    _alphaBeta( model.getAlpha() * beta ),
    _lower( _size, _size ),
    _diagonal( _size ),
    _hasBeenInitialized( false )
{}

CholeskySolver::~CholeskySolver() {
    // Deallocate memory for the Cholesky factorization
    // (no need when Blitz++ arrays are initialized as above)
}

void CholeskySolver::init() {
    // Return if CholeskySolver has already been initialized
    if ( _hasBeenInitialized ) return;

    // Build matrix M explicitly, one column at a time
    array2<double> matrixM( _size, _size );
    computeMatrixM( matrixM );

    // Compute Cholesky factorization
    computeFactorization( matrixM );
    _hasBeenInitialized = true;
}

// Compute the matrix M to be factored
void CholeskySolver::computeMatrixM( array2<double>& matrixM ) {
    BoundaryVector e( _numPoints ); // basis vector
    BoundaryVector x( _numPoints ); // x = M(e)

    cerr << "Computing the matrix for Cholesky factorization..." << flush;
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
    cerr << "done" << endl;
}

// Compute the Cholesky factorization of matrixM
//    M = L L*
// where L is lower triangular.
// Preconditions:
//      M is symmetric
// Postconditions:
//      _lower contains the strictly lower triangular part of L (no diagonal)
//      _diag  contains the diagonal elements of L
void CholeskySolver::computeFactorization( const array2<double>& matrixM ) {
    
    cerr << "Computing Cholesky factorization..." << flush;
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
    cerr << "done" << endl;
}

// Load a Cholesky decomposition from a file with name <basename>.cholesky
// Return true if successful
bool CholeskySolver::load(const string& basename) {
    string filename = basename + ".cholesky";
    cerr << "Loading Cholesky factorization from file " << filename
        << "..." << flush;
    ifstream infile( filename.c_str() );
    if ( ! infile.good() ) {
        cerr << "(failed: could not open file)" << endl;
        return false;
    }

    // read dimension of matrix in file
    int n;
    infile >> n;
    if (n != _size) {
        cerr << "(failed: wrong file size)" << endl;
        return false;
    }

    // read value of alphaBeta in file
    double alphaBeta_in;
    infile >> alphaBeta_in;
    if (alphaBeta_in != _alphaBeta) {
        cerr << "(failed: wrong timestep or Re)" << endl;
        return false;
    }
    
    // read in diagonal part
    for ( int i=0; i<_size; ++i ) {
        infile >> _diagonal(i);
    }

    // check the marker, to make sure we did not get off track
    char c;
    infile >> c;
    if (c != '#') {
        cerr << "(failed: corrupt file)" << endl;
        return false;
    }
    
    // read the lower triangular portion
    for ( int i=0; i<_size; ++i) {
        for (int j=0; j<i; ++j) {
            infile >> _lower(i,j);
        }
    }
    _hasBeenInitialized = true;
    cerr << "done" << endl;
    return true;
}

// Save a Cholesky decomposition a file with name <basename>.cholesky,
// overwriting if necessary.
// Return true if successful
// Note: saves only strictly lower triangular portion of _lower, since
// that is all that is needed for back substitution
bool CholeskySolver::save(const string& basename) {
    assert( _hasBeenInitialized );
    string filename = basename + ".cholesky";
    cerr << "Saving Cholesky factorization to file " << filename
        << "..." << flush;
    ofstream outfile( filename.c_str() );
    if ( ! outfile.good() ) {
        cerr << "(failed: could not open file)" << endl;
        return false;
    }
    outfile << _size << endl;
    outfile << setprecision(17) << _alphaBeta << endl;
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
    cerr << "done" << endl;
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
