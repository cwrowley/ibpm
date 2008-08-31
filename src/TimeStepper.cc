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

#include <string>
using namespace std;

namespace ibpm {

TimeStepper::TimeStepper(
    NavierStokesModel& model,
    double timestep) :
    _model(model),
    _timestep(timestep) {
}

ProjectionSolver* TimeStepper::createSolver(double alpha) {
    // Check whether all bodies are stationary
    //      If so, return a CholeskySolver
    //      If not, return a ConjugateGradientSolver
    
    if ( _model.getGeometry().isStationary() ) {
        cout << "Using Cholesky solver" << endl;
        return new CholeskySolver( _model, alpha );
    }
    else {
        double tol = 1e-7;
        cout << "Using ConjugateGradient solver, tolerance = " << tol << endl;
        return new ConjugateGradientSolver( _model, alpha, tol );    
    }
}

void TimeStepper::init() {
    _model.init();
}

void TimeStepper::loadState(string filename) {}

void TimeStepper::saveState(string filename) {}

} // namespace ibpm
