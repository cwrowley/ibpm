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
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "OutputForce.h"
#include "State.h"
#include "Output.h"
#include "VectorOperations.h"
#include <stdio.h>
#include <string>
using namespace std;

namespace ibpm {

OutputForce::OutputForce(string filename) :
    _filename( filename )
{}

bool OutputForce::init() {
    _fp = fopen( _filename.c_str(), "w" );
    if ( _fp == NULL ) return false;
    else return true;
}

bool OutputForce::cleanup() {
    bool status = true;
    if ( _fp != NULL ) {
        status = fclose( _fp );
    }
    return status;
}

bool OutputForce::doOutput(const State& x) {
    double drag = 0.;
    double lift = 0.;
 
    x.computeNetForce( drag, lift );

    // Convert forces to lift and drag coefficients:
    // If L_d is dimensional lift, then in the nondimensionalization of the
    // code (lengths by c, density by rho, velocity by U), we have
    //    L = L_d / (c rho U^2)
    // so
    //    C_L = L_d / (1/2 rho U^2)
    //        = 2 L
    drag *= 2;
    lift *= 2;

    if ( _fp == NULL ) return false;
    fprintf( _fp, "%5d %.5e %.5e %.5e\n", x.timestep, x.time, drag, lift );   
    fflush( _fp );
    
    return true;
}

} // namespace ibpm
