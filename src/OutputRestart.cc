// OutputRestart.cc
//
// Description:
// Implementation of output routine for writing restart files.
//
// Author(s):
// Clancy Rowley
//
// Date: 22 Aug 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "OutputRestart.h"
#include "State.h"
#include "Output.h"
#include <stdio.h>
#include <string>

namespace ibpm {

OutputRestart::OutputRestart(string formatString) {
    _formatString = formatString;
}

bool OutputRestart::doOutput(const State& x) {
    char filename[256];
    sprintf( filename, _formatString.c_str(), x.timestep );
    bool status = x.save( filename );
    return status;
}

bool OutputRestart::doOutput(const BaseFlow& q, const State& x) {
    // Currently no use for baseflow, but this method is defined for future 
    // flexibility
    return doOutput(x);
}

void OutputRestart::setFilename( string formatString ) {
    _formatString = formatString;
}

} // namespace ibpm

