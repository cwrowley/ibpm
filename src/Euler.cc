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
// $HeadURL:// $Header$

#include "BoundaryVector.h"
#include "Geometry.h"
#include "State.h"
#include "NavierStokesModel.h"
#include "ProjectionSolver.h"
#include "TimeStepper.h"
#include "Euler.h"

Euler::Euler(const NavierStokesModel& model, double timestep) :
    TimeStepper(model, timestep) {
    _solver = createSolver(timestep);
}

void Euler::advance(State& x) {
    // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
    Scalar a = x.gamma;
    
    // Evaluate Right-Hand-Side (b) for second equation of ProjecitonSolver
    const Geometry* geom = _model->getGeometry();
    BoundaryVector b = geom->getVelocities();
    
    // Call the ProjectionSolver to determine the circulation and forces
    _solver->solve( a, b, x.gamma, x.f );
    
    // Compute the corresponding flux
    _model->computeFlux( x.gamma, x.q );
}
