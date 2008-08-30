// ProjectionSolver.cc
//
// Description:
// Implementation of the ProjectionSolver class
//
// Author(s):
// Clancy Rowley
//
// Date: 3 Jul 2008
//
// $Revision: 17 $
// $LastChangedDate: 2008-07-03 01:55:25 -0400 (Thu, 03 Jul 2008) $
// $LastChangedBy: clancy $
// $HeadURL: $

#include "Scalar.h"
#include "BoundaryVector.h"
#include "NavierStokesModel.h"
#include "ProjectionSolver.h"

namespace ibpm {

// Constructor: initialize private data for eigenvalues of Ainv
ProjectionSolver::ProjectionSolver(
    const NavierStokesModel& model,
    double alpha) :
	_alpha(alpha),
	_model(&model),
    _eigenvaluesOfAinv(*( _model->getGrid() )) {

    // A = (1 - alpha/2 L)
    // Calculate eigenvalues of Ainv
    Scalar eigA = *( _model->getLambda() );
    eigA *= -_alpha/2;
    eigA += 1;
    _eigenvaluesOfAinv = 1 / eigA;
}

ProjectionSolver::~ProjectionSolver() {}

// Solve for gamma and f for a system of the form
//   A gamma + B f = a
//   C gamma = b
// using a fractional step method:
//   A gamma^* = a
//   C A^{-1}B f = C gamma^* - b
//   gamma = gamma^* - A^{-1} B f

void ProjectionSolver::solve(
    const Scalar& a,
    const BoundaryVector& b,
    Scalar& gamma,
    BoundaryVector& f
    ) {
    
    // A gamma^* = a
    Scalar gammaStar = Ainv(a);
    
    // C A^{-1}B f = C gamma^* - b
    BoundaryVector rhs = C(gammaStar);
    rhs -= b;
    Minv(rhs, f);

    // gamma = gamma^* - A^{-1} B f
    Scalar Bf = B(f);
    Scalar Ainv_B_f = Ainv( Bf );
    gamma = gammaStar - Ainv_B_f;
}

} // namespace ibpm
