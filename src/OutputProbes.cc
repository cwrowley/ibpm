// OutputProbes.cc
//
// Description:
// Implementation of output routine for writing a list of probes' states
//
// Author(s):
// Clancy Rowley
// Zhanhua Ma
//
// Date: 11 May 2008

#include "OutputProbes.h"
#include "Array.h"
#include "State.h"
#include "Output.h"
#include "VectorOperations.h"
#include <vector>
#include <stdio.h>
#include <string>
using namespace std;

namespace ibpm {

const int OutputProbes::_lev = 0;   // all probes at finest grid level
const int OutputProbes::_dimen = 2; // two-dimensional domain for now

OutputProbes::OutputProbes(string filename, Grid& grid ) :
    _filename( filename ),
    _grid( grid ),
    _info_fp( NULL )
{}

bool OutputProbes::init() {
    _info_fp = fopen( (_filename + "info").c_str(), "w" );
    if ( _info_fp == NULL ) return false;
    
    // write: probe#, probe grid indices i,  j,  probe coordinates x, y 
    for ( unsigned int n = 0 ; n < _probes.size(); n++ ) {
        int i = getProbeIndexX( n+1 );
        int j = getProbeIndexY( n+1 );
        double x = getProbeCoordX( n+1 );
        double y = getProbeCoordY( n+1 );
        fprintf( _info_fp, "%2d %3d %d %.5e %.5e\n", 
                n+1, i, j, x, y ); 
    }
    fclose( _info_fp ); 
    
    // open files for each probe to output data
    for ( unsigned int n = 0; n < _probes.size(); n++ ) {
        char name[256];
        sprintf( name, _filename.c_str(), n+1 );
        _probes[n].fp = fopen( name, "w" );
        if ( _probes[n].fp == NULL ) return false;
    }
    return true;
}

bool OutputProbes::cleanup() {
    bool status = true;
    for ( unsigned int n=0 ; n < _probes.size(); n++ ) {
        if ( _probes[n].fp != NULL ) {
                status = status && fclose( _probes[n].fp );
        }
    }
    return status;
}

bool OutputProbes::doOutput( const State& state) {
    // TODO: Unnecessary to transform velocity fields everywhere, when only a few probe points will be used
    Scalar u(_grid);
    Scalar v(_grid);
    FluxToVelocity( state.q, u, v );

    // Write  u, v, qx, qy, omega, all at gridpoint/edge (i, j), for each probe,
    // in seperate files
    // TODO: Why store qx and u?  Seems like this is redundant.
    for ( unsigned int n=0; n < _probes.size(); n++ ) {
        assert(_probes[n].fp != NULL);
        int i = _probes[n].i;
        int j = _probes[n].j;
        fprintf( _probes[n].fp, "%5d %.5e %.14e %.14e %.14e %.14e %.14e\n",  
                 state.timestep, state.time,  
                 u(_lev, i, j), v(_lev, i, j),
                 state.q(_lev, X, i, j), state.q(_lev, Y, i, j),
                 state.omega(_lev, i, j) ); 
        fflush( _probes[n].fp );
    }
    
    return true;
}

void OutputProbes::addProbeByIndex( int i, int j ){
    // check if has been initialized
    if ( _info_fp != NULL ) {
        cout << "Warning: Addition of probes is allowed only before initialization." << endl;
        exit(1);
    }
    
    if ( i < 1 || j < 1 || i > _grid.Nx()-1 || j > _grid.Ny()-1 ) {
        cout << "Warning: invalid probe position: (" << i << "," << j << ")"
            << endl
            << "Probes should be located at the inner nodes "
            << "at the finest grid level." << endl;
        exit(1);
    };
    Probe probe(i,j);
    _probes.push_back(probe);
}

void OutputProbes::addProbeByPosition( double xcord, double ycord ) {
    int i = _grid.getXGridIndex( xcord );
    int j = _grid.getYGridIndex( ycord );
    addProbeByIndex( i, j ); 
}

void OutputProbes::addProbe( int i, int j ) {
    addProbeByIndex( i, j );
}

void OutputProbes::addProbe( double xcord, double ycord ) {
    addProbeByPosition( xcord, ycord ); 
}
    
void OutputProbes::print() {
    cout << "\n-- Probe locations: by grid indices --" << endl;
    for (unsigned int n = 0; n < _probes.size(); n++ ) {
        cout << " " << n+1  << "-th probe: grid index ( " 
             << _probes[n].i << ", " << _probes[n].j << " )" << endl;
    }
    cout << "---------- end of probe list --------" << endl;
}

} // namespace ibpm
