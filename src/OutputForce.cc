// OutputForce.cc
//
// Description:
// Implementation of output routine for writing a list of force coefficients.
//
// Author(s):
// Clancy Rowley
// Steve Brunton
//
// Date: 22 Aug 2008
//
// $Revision: 117 $
// $LastChangedDate: 2008-09-01 21:43:39 -0400 (Mon, 01 Sep 2008) $
// $LastChangedBy: cwrowley $
// $HeadURL: svn+ssh://sbrunton@zion.princeton.edu/ibpm/trunk/src/OutputRestart.cc $

#include "OutputForce.h"
#include "State.h"
#include "Output.h"
#include "VectorOperations.h"
#include <stdio.h>
#include <string>
using namespace std;

namespace ibpm {

OutputForce::OutputForce(string filename) {
    _filename = filename;
}

bool OutputForce::doOutput(const State& x) {
    double drag = 0.;
    double lift = 0.;
    char filename[256];
    sprintf( filename, _filename.c_str() );
 
    computeNetForce( x.f, drag, lift);

    const Grid& grid = x.gamma.getGrid();
    drag *= 2./grid.getDx();
    lift *= 2./grid.getDx();

    FILE *fp = fopen( filename, "a" );
    if (fp == NULL) return false;
    fprintf( fp, "%0.5d %.5e %.5e %.5e\n",x.timestep, x.time, drag, lift);   
    fclose(fp);   
    return true;
}

} // namespace ibpm
