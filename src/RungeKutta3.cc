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
#include "RungeKutta3.h"
#include <stdio.h>
#include <string>
using namespace std;


namespace ibpm {

RungeKutta3::RungeKutta3(Model& model, double timestep) :
    TimeStepper("3rd-order Runge-Kutta", model, timestep),
    _linearTermEigenvalues1( model.getLambda() ),
    _linearTermEigenvalues2( model.getLambda() ),
    _linearTermEigenvalues3( model.getLambda() ), 
    _x1( model.getGrid(), model.getGeometry().getNumPoints() ),
    _x2( model.getGrid(), model.getGeometry().getNumPoints() ),
    _Q1( model.getGrid() ),
    _Q2( model.getGrid() ),
    _Q3( model.getGrid() ),
    _a( model.getGrid() ),
    _b( model.getGeometry().getNumPoints() ),
    _b0( model.getGeometry().getNumPoints() )
    {
    // Set values of constants given in Peyret, p.149[3]  
    _A1 = 0.; 
    _A2 = -5./9.;
    _A3 = -153./128.;
    _B1 = 1./3.;
    _B2 = 15./16.;
    _B3 = 8./15.;
    _Bp1 = 1./6.;
    _Bp2 = 5./24.;
    _Bp3 = 1./8.;
    // compute eigenvalues of linear terms on RHS of each projection solve
    _linearTermEigenvalues3 *= timestep * _Bp3;
    _linearTermEigenvalues3 += 1; 
    _linearTermEigenvalues2 *= timestep * _Bp2;
    _linearTermEigenvalues2 += 1;
    _linearTermEigenvalues1 *= timestep * _Bp1;
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
        _solver1 = createSolver(_timestep * 2. * _Bp1);
    }

    // 2nd solver: alpha = 5*h/12
    if (_solver2 == 0) {
        _solver2 = createSolver(_timestep * 2. * _Bp2);
    }

    // 3rd solver: alpha = h/4
    if (_solver3 == 0) {
        _solver3 = createSolver(_timestep * 2. * _Bp3);
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

    //*****
    // First Projection Solve (intermediate state _x1) 

    // RHS for 1st eqn of ProjectionSolver
    _Q1 = _model.N( x );
    _Q1 *= _timestep;

    _a = _model.S( x.gamma );
    _a *= _linearTermEigenvalues1;
    _a = _model.Sinv( _a );

    _a += _Q1 * _B1;

    // RHS for 2nd eqn of ProjectionSolver
    _b = geom.getVelocities();
    _b0 = _model.getBaseFlowBoundaryVelocities();
    _b -= _b0;
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver1->solve( _a, _b, _x1.gamma, _x1.f );
    
    // Compute the corresponding flux
    _model.computeFlux( _x1.gamma, _x1.q );

    //*****
    // Second Projection Solve (intermediate state _x2)

    // RHS for 1st eqn of ProjectionSolver
    _Q2 = _model.N( _x1 );
    _Q2 *= _timestep;
    _Q2 += _Q1 * _A2;

    _a = _model.S( _x1.gamma );
    _a *= _linearTermEigenvalues2;  
    _a = _model.Sinv( _a );

    _a += _Q2 * _B2;
    
    // RHS for 2nd eqn of ProjectionSolver
    _b = geom.getVelocities();
    _b0 = _model.getBaseFlowBoundaryVelocities();
    _b -= _b0;
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver2->solve( _a, _b, _x2.gamma, _x2.f );
    
    // Compute the corresponding flux
    _model.computeFlux( _x2.gamma, _x2.q );

    //*****
    // Third Projection Solve (final state)

    // RHS for 1st eqn of ProjectionSolver  
    _Q3 = _model.N( _x2 );
    _Q3 *= _timestep;
    _Q3 += _Q2 * _A3;

    _a = _model.S( _x2.gamma );
    _a *= _linearTermEigenvalues3;  
    _a = _model.Sinv( _a );

    _a += _Q3 * _B3;
    
    // RHS for 2nd eqn of ProjectionSolver
    _b = geom.getVelocities();
    _b0 = _model.getBaseFlowBoundaryVelocities();
    _b -= _b0;
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver3->solve( _a, _b, x.gamma, x.f );
    
    // Compute the corresponding flux
    _model.computeFlux( x.gamma, x.q );

    // Update the value of the time
    x.time += _timestep;
    ++x.timestep;
}

} // namespace ibpm
