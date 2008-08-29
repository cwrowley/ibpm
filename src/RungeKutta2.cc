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
// $Revision: 105 $
// $LastChangedDate: 2008-08-27 23:14:34 -0400 (Thu, 28 Aug 2008) $
// $LastChangedBy: sbrunton $
// $HeadURL:// $Header$

#include "BoundaryVector.h"
#include "Geometry.h"
#include "Grid.h"
#include "NavierStokesModel.h"
#include "ProjectionSolver.h"
#include "TimeStepper.h"
#include "State.h"
#include "RungeKutta2.h"

RungeKutta2::RungeKutta2(const NavierStokesModel& model, double timestep) :
    TimeStepper(model, timestep),
    _linearTermEigenvalues( *(_model->getLambda()) ) {
    // compute eigenvalues of linear term on RHS
    _linearTermEigenvalues *= timestep/2.;
    _linearTermEigenvalues += 1;
  
    _grid = _model->getGrid(); 
    _geom = _model->getGeometry();
    _solver = createSolver(timestep);
}

void RungeKutta2::advance(State& x) {
    State x1( *(_grid), *(_geom));

    // If the body is moving, update the positions of the bodies
    if ( ! _geom->isStationary() ) {
        _geom->moveBodies(x.time);
    }

    // First Projection Solve 

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
    Scalar a = _model->S( x.gamma );
    a *= _linearTermEigenvalues;
    a = _model->Sinv( a );

    Scalar nonlinear = _model->nonlinear( x );
    a += _timestep * nonlinear;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    BoundaryVector b = _geom->getVelocities();
    BoundaryVector b0 = _model->getBaseFlowBoundaryVelocities();
    b -= b0;
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver->solve( a, b, x1.gamma, x1.f );
    
    // Compute the corresponding flux
    _model->computeFlux( x1.gamma, x1.q );

    // Second Projection Solve

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
    a = _model->S( x.gamma );
    a *= _linearTermEigenvalues;
    a = _model->Sinv( a );

    nonlinear = _model->nonlinear( x ) + _model->nonlinear( x1 );
    a += .5 * _timestep * nonlinear;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    b = _geom->getVelocities();
    b0 = _model->getBaseFlowBoundaryVelocities();
    b -= b0;
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver->solve( a, b, x.gamma, x.f );
    
    // Compute the corresponding flux
    _model->computeFlux( x.gamma, x.q );
 
   
    // Update the value of the time
    x.time += _timestep;
    ++x.timestep;
}
