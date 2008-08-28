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

// Allocate memory for the Cholesky factorization, but do not compute it
CholeskySolver::CholeskySolver(
    const NavierStokesModel& model,
    double alpha
    ) :
    ProjectionSolver( model, alpha ) {
    // Determine size of the system to be solved
    // Allocate memory for the Cholesky factorization
}

CholeskySolver::~CholeskySolver() {
    // Deallocate memory for the Cholesky factorization
}

// Compute the Cholesky factorization
void CholeskySolver::computeCholesky() {
    // Build matrix M explicitly, one column at a time
    // Compute Cholesky factorization
}

// Load a Cholesky decomposition from the specified file.
// Return true if successful
bool CholeskySolver::loadCholesky(const string& filename) {
    return false;
}

// Save a Cholesky decomposition to the specified file, overwriting if
// necessary.
// Return true if successful
bool CholeskySolver::saveCholesky(const string& filename) {
    return false;
}

// Solve A x = b using the Cholesky factorization A = LL*
void CholeskySolver::Minv(
    const BoundaryVector& b,
    BoundaryVector& x
    ) {
    // TODO: Work here
    // Solve L y = b for y
    // Solve L^Tx = y for x
    
}
