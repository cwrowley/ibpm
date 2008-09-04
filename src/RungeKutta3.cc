// RungeKutta3.cc
//
// Description:
// Implementation of the RungeKutta3 class
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
#include "RungeKutta3.h"
#include <stdio.h>
#include <string>
using namespace std;


namespace ibpm {

RungeKutta3::RungeKutta3(NavierStokesModel& model, double timestep) :
    TimeStepper(model, timestep),
    _linearTermEigenvalues1( model.getLambda() ),
    _linearTermEigenvalues2( model.getLambda() ),
    _linearTermEigenvalues3( model.getLambda() ), 
    _x1( model.getGrid(), model.getGeometry() ),
    _x2( model.getGrid(), model.getGeometry() )
    {
    // compute eigenvalues of linear terms on RHS of each projection solve
    _linearTermEigenvalues3 *= timestep/8.; 
    _linearTermEigenvalues3 += 1; 
    _linearTermEigenvalues2 *= timestep*5./24.;
    _linearTermEigenvalues2 += 1;
    _linearTermEigenvalues1 *= timestep/6.;
    _linearTermEigenvalues1 += 1;
    _solver1 = 0; // initialize solvers in init()
    _solver2 = 0;
    _solver3 = 0;
}

RungeKutta3::~RungeKutta3() {
    delete _solver1;
    delete _solver2;
    delete _solver3;
}

void RungeKutta3::createAllSolvers() {
    // 1st solver: alpha = h/3
    if (_solver1 == 0) {
        _solver1 = createSolver(_timestep/3.);
    }

    // 2nd solver: alpha = 5*h/12
    if (_solver2 == 0) {
        _solver2 = createSolver(_timestep*5./12.);
    }

    // 3rd solver: alpha = h/4
    if (_solver3 == 0) {
        _solver3 = createSolver(_timestep/4.);
    }
}

void RungeKutta3::init() {
    createAllSolvers();
    _solver1->init();
    _solver2->init();
    _solver3->init();
}
 
bool RungeKutta3::load(const string& basename) {
    createAllSolvers();
    bool success;
    success =            _solver1->load( basename + "_01" );
    success = success && _solver2->load( basename + "_02" );
    success = success && _solver3->load( basename + "_03" );
    return success;
}

// Save the ProjectionSolver's state, if necessary
// (e.g. the Cholesky factorization)
bool RungeKutta3::save(const string& basename) {
    bool success;
    // save first solver
    if (_solver1 != 0) {
        success = _solver1->save( basename + "_01" );
    }
    // save second solver
    if (_solver2 != 0) {
        success = success && _solver2->save( basename + "_02" );
    }
    // save third solver
    if (_solver3 != 0) {
        success = success && _solver3->save( basename + "_03" );
    }
    return success;
}

void RungeKutta3::advance(State& x) {
    // If the body is moving, update the positions of the bodies
    const Geometry& geom = _model.getGeometry();
    if ( ! geom.isStationary() ) {
        geom.moveBodies(x.time);
    }

    // First Projection Solve 

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
    Scalar Q1 = _model.nonlinear( x );
    Q1 *= _timestep;

    Scalar a = _model.S( x.gamma );
    a *= _linearTermEigenvalues1;
    a = _model.Sinv( a );

    a += Q1 / 3.;

    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    BoundaryVector b = geom.getVelocities();
    BoundaryVector b0 = _model.getBaseFlowBoundaryVelocities();
    b -= b0;
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver1->solve( a, b, _x1.gamma, _x1.f );
    
    // Compute the corresponding flux
    _model.computeFlux( _x1.gamma, _x1.q );

    // Second Projection Solve

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver  
    Scalar Q2 = _model.nonlinear( _x1 );
    Q2 *= _timestep;
    Q2 += -Q1 * 5. / 9.;

    //a = _model.S( _x1.gamma );
    a = _model.S( _x1.gamma );
    a *= _linearTermEigenvalues2;  
    a = _model.Sinv( a );

    a += Q2 * 15. / 16.;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    b = geom.getVelocities();
    b0 = _model.getBaseFlowBoundaryVelocities();
    b -= b0;
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver2->solve( a, b, _x2.gamma, _x2.f );
    
    // Compute the corresponding flux
    _model.computeFlux( _x2.gamma, _x2.q );

    // Third Projection Solve

    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver  
    Scalar Q3 = _model.nonlinear( _x2 );
    Q3 *= _timestep;
    Q3 += -Q2 * 153. / 128.;

    //a = _model.S( _x2.gamma );
    a = _model.S( _x2.gamma );
    a *= _linearTermEigenvalues3;  
    a = _model.Sinv( a );

    a += Q3 * 8. / 15.;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
    b = geom.getVelocities();
    b0 = _model.getBaseFlowBoundaryVelocities();
    b -= b0;
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver2->solve( a, b, x.gamma, x.f );
    
    // Compute the corresponding flux
    _model.computeFlux( x.gamma, x.q );

    // Update the value of the time
    x.time += _timestep;
    ++x.timestep;
}

} // namespace ibpm
