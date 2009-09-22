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
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "Scalar.h"
#include "BoundaryVector.h"
#include "NavierStokesModel.h"
#include "ProjectionSolver.h"
#include <string>

namespace ibpm {

// Constructor: initialize Helmholtz solver to solve (1 - beta/2 L) u = f
// 
ProjectionSolver::ProjectionSolver(
    const Grid& grid,
    const NavierStokesModel& model,
    double beta) :
    _beta(beta),
    _grid(grid),
    _model(model),
    _helmholtz( grid, -beta/2. * _model.getAlpha() )
{}

ProjectionSolver::~ProjectionSolver() {}

// Initialization for this base class done in constructor
// Subclasses use this for their own initialization, if needed
void ProjectionSolver::init() {}

// Subclasses can override these to save and load their own state
// (e.g. a Cholesky factorization)
// By default, not implemented: return false
bool ProjectionSolver::save(const std::string& filename) { return false; }
bool ProjectionSolver::load(const std::string& filename) { return false; }


// Solve for omega and f for a system of the form
//   A omega + B f = a
//   C omega = b
// using a fractional step method:
//   A omega^* = a
//   C A^{-1}B f = C omega^* - b
//   omega = omega^* - A^{-1} B f

void ProjectionSolver::solve(
    const Scalar& a,
    const BoundaryVector& b,
    Scalar& omega,
    BoundaryVector& f
    ) {
    
    // A omega^* = a
    Scalar omegaStar( a.getGrid() );
    Ainv( a, omegaStar );
    
    // C A^{-1}B f = C omega^* - b
    BoundaryVector rhs( f.getNumPoints() );
    C( omegaStar, rhs );
    rhs -= b;               // rhs = C omega^* - b
    Minv( rhs, f );         // f = Minv( rhs )

    // omega = omega^* - A^{-1} B f
    Scalar c( a.getGrid() );
    B( f, c );              // c = Bf
    Ainv( c, c );           // c = Ainv(Bf)
    omega = omegaStar - c;
}
    
void ProjectionSolver::Ainv(const Scalar& x, Scalar& y) {
    _helmholtz.solve( x, y );
}
    
void ProjectionSolver::B(const BoundaryVector& f, Scalar& y) {
    _model.B( f, y );
    y *= _beta;
}


void ProjectionSolver::C(const Scalar& x, BoundaryVector& y ) {
    _model.C( x, y );
}

// Compute y = M(f), where M = C A^{-1} B
void ProjectionSolver::M(const BoundaryVector& f, BoundaryVector& y ) {
    Scalar u( _grid );
    B( f, u );          // u = B f
    Ainv( u, u );       // u = Ainv B f
    C( u, y );          // y = C Ainv B f
}    

BoundaryVector ProjectionSolver::M(const BoundaryVector& f) {
    BoundaryVector y( f.getNumPoints() );
    M( f, y );
    return y;
}
    
} // namespace ibpm
