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
// $Revision:  $
// $LastChangedDate:  $
// $LastChangedBy: clancy $
// $HeadURL: $

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
    char buf[256];
    sprintf( buf, _formatString.c_str(), x.timestep );

    // TODO: Work here
    // Open the file
    // Write contents to file
    // Close the file
    cout << "Writing restart file " << buf << endl;
    cout << " NOTE: Not actually implemented!" << endl;
    return true;
}

} // namespace ibpm
