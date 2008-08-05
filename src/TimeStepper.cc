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
// $HeadURL:// $Header$

#include "NavierStokesModel.h"
#include "Geometry.h"

TimeStepper(const NavierStokesModel& model, double timestep) {}

ProjectionSolver* TimeStepper::createSolver(double alpha) {
    // TODO: Check whether all bodies are stationary
    //      If so, return a CholeskySolver
    //      If not, return a ConjugateGradientSolver
    
    // For now, just return a ConjugateGradientSolver
    return new ConjugateGradientSolver(alpha);
}
