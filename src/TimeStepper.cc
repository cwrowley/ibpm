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

#include "Model.h"
#include "Geometry.h"
#include "ProjectionSolver.h"
#include "ConjugateGradientSolver.h"
#include "CholeskySolver.h"
#include "TimeStepper.h"

#include <string>
using namespace std;

namespace ibpm {

TimeStepper::TimeStepper(
    string name,
    Grid& grid,
    Model& model,
    double timestep) :
    _name(name),
    _grid(grid),
    _model(model),
    _timestep(timestep) {
}

string TimeStepper::getName() {
    return _name;
}

ProjectionSolver* TimeStepper::createSolver(double beta) {
    // Check whether all bodies are stationary
    //      If so, return a CholeskySolver
    //      If not, return a ConjugateGradientSolver
    
    if ( _model.isTimeDependent() ) {
        double tol = 1e-7;
        cerr << "Using ConjugateGradient solver for projection step" << endl
        << "  tolerance = " << tol << endl;
        return new ConjugateGradientSolver( _grid, _model, beta, tol );    
    }
    else {
        cerr << "Using Cholesky solver for projection step" << endl;
        return new CholeskySolver( _grid, _model, beta );
    }
}

void TimeStepper::init() {
    _model.init();
}

// Implemented in subclasses: return false by default
bool TimeStepper::load(const string& filename) { return false; }

// Implemented in subclasses: return false by default
bool TimeStepper::save(const string& filename) { return false; }

} // namespace ibpm
