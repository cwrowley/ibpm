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
using namespace std;

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

void OutputRestart::setFilename( string formatString ) {
    _formatString = formatString;
}

} // namespace ibpm

