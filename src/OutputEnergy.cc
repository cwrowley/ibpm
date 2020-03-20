#include "OutputEnergy.h"
#include "State.h"
#include "Output.h"
#include "VectorOperations.h"
#include <stdio.h>
#include <string>

namespace ibpm {

OutputEnergy::OutputEnergy(string filename) :
    _filename( filename )
{}

bool OutputEnergy::init() {
    _fp = fopen( _filename.c_str(), "w" );
    if ( _fp == NULL ) return false;
    else return true;
}

bool OutputEnergy::cleanup() {
    bool status = true;
    if ( _fp != NULL ) {
        status = fclose( _fp );
    }
    return status;
}
    
bool OutputEnergy::doOutput(const State& x) {
    double energy = 0.;
    energy = .5 * InnerProduct( x.q, x.q );
    
    if ( _fp == NULL ) return false;
    fprintf( _fp, "%5d %.5e %.5e\n", x.timestep, x.time, energy );   
    fflush( _fp );
    return true;
}
    
bool OutputEnergy::doOutput(const BaseFlow & q, const State& x) {
    // Currently no use for baseflow, but this method is defined for future 
    // flexibility
    return doOutput(x);
}
    
} // namespace ibpm
