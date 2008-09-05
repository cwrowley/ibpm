// AdamsBashforth.cc
//
// Description:
// Implementation of the AdamsBashforth class
//
// Author(s):
// Clancy Rowley
// Steve Brunton
//
// Date: 28 Aug 2008
//
// $Revision: 105 $
// $LastChangedDate: 2008-08-27 23:14:34 -0400 (Thu, 28 Aug 2008) $
// $LastChangedBy: sbrunton $
// $HeadURL: Header $

#include "BoundaryVector.h"
#include "Geometry.h"
#include "Grid.h"
#include "NavierStokesModel.h"
#include "ProjectionSolver.h"
#include "TimeStepper.h"
#include "State.h"
#include "AdamsBashforth.h"

namespace ibpm {

AdamsBashforth::AdamsBashforth(NavierStokesModel& model, double timestep) :
    TimeStepper("Adams Bashforth", model, timestep),
    _linearTermEigenvalues( model.getLambda() ),
    _xold( model.getGrid(), model.getGeometry() ),
    _xtemp( model.getGrid(), model.getGeometry() )
    {
    // compute eigenvalues of linear term on RHS
    _linearTermEigenvalues *= timestep/2.;
    _linearTermEigenvalues += 1;
    _solver = 0; // initialize solver in init()
    _oldSaved = false;
}

AdamsBashforth::~AdamsBashforth() {
    delete _solver;
}

void AdamsBashforth::init() {
    if (_solver == 0) {
        _solver = createSolver(_timestep);
    }
    _solver->init();
}

bool AdamsBashforth::load(const string& basename) {
    if (_solver == 0) {
        _solver = createSolver(_timestep);
    }
    bool success = _solver->load( basename );
    return success;
}

// Save the ProjectionSolver's state, if necessary
// (e.g. the Cholesky factorization)
bool AdamsBashforth::save(const string& basename) {
    bool success = false;
    if (_solver != 0) {
        success = _solver->save( basename );
    }
    return success;
}

void AdamsBashforth::advance(State& x) {
    // If the body is moving, update the positions of the bodies
    const Geometry& geom = _model.getGeometry();
    if ( ! geom.isStationary() ) {
        geom.moveBodies(x.time);
    }

    // Keep unmodified state to save as previous at the end of advance 
    setTempState( x );  

    // Initialize _xold with current state if not previously saved
    if ( _oldSaved == false ) {
        setPreviousState( x );    
    }

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
    Scalar a = _model.S( x.gamma );
    a *= _linearTermEigenvalues;
    a = _model.Sinv( a );

    Scalar nonlinear = _model.nonlinear( x );
    nonlinear *= 3.;
    nonlinear -= _model.nonlinear( _xold );
    a += .5 * _timestep * nonlinear;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    BoundaryVector b = geom.getVelocities();
    BoundaryVector b0 = _model.getBaseFlowBoundaryVelocities();
    b -= b0;
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver->solve( a, b, x.gamma, x.f );
    
    // Compute the corresponding flux
    _model.computeFlux( x.gamma, x.q );

    setPreviousState( _xtemp ); 

    // Update the value of the time
    x.time += _timestep;
    ++x.timestep;
}

} // namespace ibpm
