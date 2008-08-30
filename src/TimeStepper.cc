// TimeStepper.cc
//
// Description:
// Implementation of the TimeStepper abstract base class
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

#include "NavierStokesModel.h"
#include "Geometry.h"
#include "ProjectionSolver.h"
#include "ConjugateGradientSolver.h"
#include "CholeskySolver.h"
#include "TimeStepper.h"

namespace ibpm {

TimeStepper::TimeStepper(
    const NavierStokesModel& model,
    double timestep) :
    _model(&model),
    _timestep(timestep) {
}

ProjectionSolver* TimeStepper::createSolver(double alpha) {
    // Check whether all bodies are stationary
    //      If so, return a CholeskySolver
    //      If not, return a ConjugateGradientSolver
    
    if ( _model->getGeometry()->isStationary() ) {
        return new CholeskySolver( *_model, alpha );
    }
    else {
        double tol = 1e-7;
        return new ConjugateGradientSolver( *_model, alpha, tol );    
    }
}

} // namespace ibpm
