#include "OutputEnergy.h"
#include "State.h"
#include "Output.h"
#include "VectorOperations.h"
#include <stdio.h>
#include <string>
using namespace std;

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
    const Grid& grid = x.omega.getGrid();
    Scalar u(grid);
    Scalar v(grid);
    FluxToXVelocity(x.q,u);
    FluxToYVelocity(x.q,v);
    energy = .5 * sqrt( InnerProduct( u, u ) + InnerProduct( v, v ) );

    if ( _fp == NULL ) return false;
    fprintf( _fp, "%5d %.5e %.5e\n", x.timestep, x.time, energy );   
    fflush( _fp );
    return true;
}

} // namespace ibpm
