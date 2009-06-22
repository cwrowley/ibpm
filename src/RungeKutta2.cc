// RungeKutta2.cc
//
// Description:
// Implementation of the RungeKutta2 class
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
#include "RungeKutta2.h"

namespace ibpm {

RungeKutta2::RungeKutta2( Grid& grid, Model& model, double timestep ) :
    TimeStepper("2nd-order Runge-Kutta", grid, model, timestep),
    _x1( grid, model.getNumPoints() )
    {
    _solver = createSolver( _timestep );
    // initialize solver in init()
}

RungeKutta2::~RungeKutta2() {
    delete _solver;
}

void RungeKutta2::init() {
    _solver->init();
}

bool RungeKutta2::load(const string& basename) {
    bool success = _solver->load( basename );
    return success;
}

// Save the ProjectionSolver's state, if necessary
// (e.g. the Cholesky factorization)
bool RungeKutta2::save(const string& basename) {
    bool success = false;
    if (_solver != 0) {
        success = _solver->save( basename );
    }
    return success;
}

void RungeKutta2::advance(State& x) {
    // Update the value of the time
    x.time += _timestep;
    ++x.timestep;
	
    // If the body is moving, update the positions of the bodies
    if ( _model.isTimeDependent() ) {
        _model.updateOperators( x.time );
    }

    // First Projection Solve 

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
    Scalar a = Laplacian( x.omega );
    a *= 0.5 * _timestep * _model.getAlpha();
    a += x.omega;

    Scalar nonlinear = _model.N( x );
    a += _timestep * nonlinear;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    BoundaryVector b = _model.getConstraints();
    
    // Call the ProjectionSolver to determine the vorticity and forces
    _solver->solve( a, b, _x1.omega, _x1.f );
    
    // Update the rest of the state (i.e. flux)
    _model.refreshState( _x1 );

    // Second Projection Solve

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
    // TODO: Redo this to avoid calculating L(omega) and N(x) twice 
    a = Laplacian( x.omega );
    a *= 0.5 * _timestep * _model.getAlpha();
    a += x.omega;

    nonlinear = _model.N( x ) + _model.N( _x1 );
    a += 0.5 * _timestep * nonlinear;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    b = _model.getConstraints();
    
    // Call the ProjectionSolver to determine the vorticity and forces
    _solver->solve( a, b, x.omega, x.f );
    
    // Update the rest of the state (i.e. flux)
    _model.refreshState( x );

}

} // namespace ibpm
