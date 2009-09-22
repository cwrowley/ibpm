#ifndef _IBPM_H_
#define _IBPM_H_

/*!
    \file ibpm.h
    
    \brief Main header file for immersed boundary projection method (IBPM) library
    
    \author Clancy Rowley
    \author $LastChangedBy$
    \date 30 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

#define IBPM_VERSION "1.0"

// basic classes
#include "Grid.h"
#include "RigidBody.h"
#include "Geometry.h"

// data structures
#include "Scalar.h"
#include "Flux.h"
#include "BoundaryVector.h"
#include "State.h"

// operations
#include "VectorOperations.h"
#include "NavierStokesModel.h"

// timesteppers
#include "IBSolver.h"

// motion
#include "Motion.h"
#include "FixedPosition.h"
#include "PitchPlunge.h"

// output routines
#include "Logger.h"
#include "OutputTecplot.h"
#include "OutputRestart.h"
#include "OutputEnergy.h"
#include "OutputForce.h"
#include "OutputProbes.h"

// utilities
#include "utils.h"
#include "ParmParser.h"

#endif /* _IBPM_H_ */
