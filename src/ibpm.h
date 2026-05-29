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
#include "BaseFlow.h"
#include "State.h"
#include "StateVector.h"

// operations
#include "VectorOperations.h"
#include "NavierStokesModel.h"

// timesteppers
#include "IBSolver.h"

// motion
#include "Motion.h"
#include "FixedPosition.h"
#include "FixedVelocity.h"
#include "PitchPlunge.h"
// #include "PitchPlungeT.h"
// #include "PitchPlungeFS.h"
// #include "PitchPlungeFS_test.h"
// #include "PitchPlungeT_JF.h"
// #include "PitchPlungeT_adj_check.h"
// #include "PitchPlungeT_adj_check_full.h"
#include "PitchPlungeT_mod.h"
#include "SigmoidalStep.h"
#include "LagStep1.h"
#include "LagStep2.h"
#include "EldredgeManeuver.h"
#include "EldredgeCombined2.h" // combines EldredgeManeuver for pitch and Eldredge2 for plunge
#include "Eldredge1.h"
#include "Eldredge2.h"
#include "MotionFile.h"
#include "MotionFilePeriodic.h"

// output routines
#include "Logger.h"
#include "OutputTecplot.h"
#include "OutputRestart.h"
#include "OutputEnergy.h"
#include "OutputForce.h"
#include "OutputProbes.h"
#include "ScalarToTecplot.h"

// utilities
#include "utils.h"
#include "ParmParser.h"

#endif /* _IBPM_H_ */
