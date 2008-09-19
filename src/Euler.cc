// Euler.cc
//
// Description:
// Implementation of the Euler class
//
// Author(s):
// Clancy Rowley
//
// Date: 2 Aug 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "BoundaryVector.h"
#include "Geometry.h"
#include "State.h"
#include "Model.h"
#include "ProjectionSolver.h"
#include "TimeStepper.h"
#include "VectorOperations.h"
#include "Euler.h"

namespace ibpm {

Euler::Euler(Grid& grid, Model& model, double timestep) :
    TimeStepper("Explicit Euler", grid, model, timestep)
    {
    _solver = createSolver( _timestep );
    // initialize solver in init()
}

Euler::~Euler() {
    delete _solver;
}

void Euler::init() {
    _solver->init();
}

bool Euler::load(const string& basename) {
    bool success = _solver->load( basename );
    return success;
}

// Save the ProjectionSolver's state, if necessary
// (e.g. the Cholesky factorization)
bool Euler::save(const string& basename) {
    bool success = false;
    if (_solver != 0) {
        success = _solver->save( basename );
    }
    return success;
}

void Euler::advance(State& x) {
    // Update the value of the time
    x.time += _timestep;
    ++x.timestep;

    // If the body is moving, update the positions of the bodies
    if ( _model.isTimeDependent() ) {
        _model.updateOperators( x.time );
    }

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
    Scalar a = Laplacian( x.gamma );
    a *= 0.5 * _timestep * _model.getAlpha();
    a += x.gamma;

    Scalar nonlinear = _model.N( x );
    a += _timestep * nonlinear;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    BoundaryVector b = _model.getConstraints();
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver->solve( a, b, x.gamma, x.f );
    
    // Update the state, for instance to compute the corresponding flux
    _model.refreshState( x );

}

} // namespace ibpm
