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
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "BoundaryVector.h"
#include "Geometry.h"
#include "Grid.h"
#include "Model.h"
#include "ProjectionSolver.h"
#include "TimeStepper.h"
#include "State.h"
#include "VectorOperations.h"
#include "AdamsBashforth.h"

namespace ibpm {

AdamsBashforth::AdamsBashforth(Grid& grid, Model& model, double timestep) :
    TimeStepper("Adams Bashforth", grid, model, timestep),
    _xold( grid, model.getNumPoints() ),
    _xtemp( grid, model.getNumPoints() )
    {
    _solver = createSolver( _timestep);
    // initialize solver in init()
    _oldSaved = false;
}

AdamsBashforth::~AdamsBashforth() {
    delete _solver;
}

void AdamsBashforth::init() {
    _solver->init();
}

bool AdamsBashforth::load(const string& basename) {
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
    // Update the value of the time
    x.time += _timestep;
    ++x.timestep;

    // If the body is moving, update the positions of the bodies
    if ( _model.isTimeDependent() ) {
        _model.updateOperators(x.time);
    }

    // Keep unmodified state to save as previous at the end of advance 
    setTempState( x );  

    // Initialize _xold with current state if not previously saved
    if ( _oldSaved == false ) {
        setPreviousState( x );    
    }

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
    Scalar a = Laplacian( x.gamma );
    a *= 0.5 * _timestep * _model.getAlpha();
    a += x.gamma;

    Scalar nonlinear = _model.N( x );
    nonlinear *= 3.;
    nonlinear -= _model.N( _xold );
    a += .5 * _timestep * nonlinear;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    BoundaryVector b = _model.getConstraints();
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver->solve( a, b, x.gamma, x.f );
    
    // Compute the corresponding flux
    _model.refreshState( x );

    setPreviousState( _xtemp ); 
}

} // namespace ibpm
