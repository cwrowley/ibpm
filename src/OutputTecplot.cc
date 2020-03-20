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

namespace ibpm {

OutputTecplot::OutputTecplot( string filename, string title, bool TecplotAllGrids ) {
    _filename = filename;
    _title = title;
    _TecplotAllGrids = TecplotAllGrids;
}
    
bool OutputTecplot::doOutput(const State& state) {
    // Add timestep to filename and title
    char filename[256];
    sprintf( filename, _filename.c_str(), state.timestep );
    char title[256];
    sprintf( title, _title.c_str(), state.timestep );
    bool status = false;    
    const Grid& grid = state.omega.getGrid();	
		
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
    if(_TecplotAllGrids) {
        status = true;
        for(int i=0;i < grid.Ngrid(); i++) {
            sprintf( filename, _filename.c_str(), state.timestep, i );
            cout << filename << endl;
            status = status && ScalarToTecplot( varVec, varNameVec, filename, title, i);
        }
    }
    else status = ScalarToTecplot( varVec, varNameVec, filename, title );
    return status;
}
    
bool OutputTecplot::doOutput(const BaseFlow& q, const State& x) {
    // Currently no use for baseflow, but this method is defined for future 
    // flexibility
    return doOutput(x);
}
    
void OutputTecplot::setFilename( string filename ) {
    _filename = filename;
}
    
void OutputTecplot::setTitle( string title ) {
    _title = title;
}
    
} // namespace ibpm
