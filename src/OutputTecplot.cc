// OutputTecplot.cc
//
// Description:
// Implementation of output routine for writing Tecplot ASCII files.
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

#include "OutputTecplot.h"
#include "State.h"
#include "Output.h"
#include "VectorOperations.h"
#include "ScalarToTecplot.h"
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

namespace ibpm {

OutputTecplot::OutputTecplot( string filename, string title ) {
    _filename = filename;
    _title = title;
}
    
bool OutputTecplot::doOutput(const BaseFlow& q, const State& state) {
    // Add timestep to filename and title
    char filename[256];
    sprintf( filename, _filename.c_str(), state.timestep );
    char title[256];
    sprintf( title, _title.c_str(), state.timestep );
    			
    // Calculate velocities
    Scalar u( state.omega.getGrid() );
    Scalar v( state.omega.getGrid() );
    FluxToVelocity( state.q, u, v );
    
    // Create vector of Scalar fields
    vector<const Scalar*> varVec;
    varVec.push_back( &u );
    varVec.push_back( &v);
    varVec.push_back( &state.omega );
    
    vector<string> varNameVec;
    varNameVec.push_back( "u" );
    varNameVec.push_back( "v" );
    varNameVec.push_back( "Vorticity" );
        
    // Write the tecplot file
    bool status = ScalarToTecplot( varVec, varNameVec, filename, title );
    return status;
}
    
void OutputTecplot::setFilename( string filename ) {
    _filename = filename;
}
    
void OutputTecplot::setTitle( string title ) {
    _title = title;
}
    
} // namespace ibpm
