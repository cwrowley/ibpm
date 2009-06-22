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
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

namespace ibpm {

OutputTecplot::OutputTecplot( string filename, string title ) {
    _filename = filename;
    _title = title;
}

// Small class for storing a list of pointers to variables, and their names
class VarList {
public:
    void addVariable( const Scalar* var, string name ) {
        _vars.push_back( var );
        _names.push_back( name );
    }

    int getNumVars() const {
        return _vars.size();
    }

    string getName(int i) const {
        return _names[i];
    }
    
    const Scalar* getVariable(int i) const {
        return _vars[i];
    }
    
private:
    vector<const Scalar*> _vars;
    vector<string> _names;
};

// Opens a new file with the specified name, and writes an ASCII Tecplot file
// containing the variables in the VarList passed in
bool writeTecplotFileASCII(
    const char* filename,
    const char* title,
    const VarList& list
);

bool OutputTecplot::doOutput(const State& state) {
    // Add timestep to filename and title
    char filename[256];
    sprintf( filename, _filename.c_str(), state.timestep );
    char title[256];
    sprintf( title, _title.c_str(), state.timestep );
    			
    // Get grid dimensions
    const Grid& grid = state.omega.getGrid();
    int nx = grid.Nx();
    int ny = grid.Ny();
    int ngrid = grid.Ngrid();

    // Calculate the variables for output
    // Calculate the grid
    Scalar x(grid);
    Scalar y(grid);
    for (int lev=0; lev<ngrid; ++lev) {
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny; ++j) {
                x(lev,i,j) = grid.getXEdge(lev,i);
                y(lev,i,j) = grid.getYEdge(lev,j);
            }
        }
    }
    
    // Calculate velocities
    Scalar u(grid);
    Scalar v(grid);
    FluxToVelocity( state.q, u, v );
        
    // Store pointers to variables and corresponding names in vectors
    VarList list;
    list.addVariable( &x, "x" );
    list.addVariable( &y, "y" );
    list.addVariable( &u, "u" );
    list.addVariable( &v, "v" );
    list.addVariable( &(state.omega), "Vorticity" );

    // Write the Tecplot file
    bool status = writeTecplotFileASCII( filename, title, list );
    return status;
}

// Write a Tecplot file with the specified filename and title,
// with data given in the VarList passed in
// For now, write only finest grid (level 0)
bool writeTecplotFileASCII(
    const char* filename,
    const char* title,
    const VarList& list
    ) {

    int numVars = list.getNumVars();
    assert( numVars > 0 );
    // Get grid information
    const Grid& grid = list.getVariable(0)->getGrid();
    int nx = grid.Nx();
    int ny = grid.Ny();

    // Write the header for the Tecplot file
    cerr << "Writing Tecplot file " << filename << endl;
    FILE *fp = fopen( filename, "w" );
    if (fp == NULL) return false;
    fprintf( fp, "TITLE = \"%s\"\n", title );
    fprintf( fp, "VARIABLES = ");
    for (int i=0; i<numVars; ++i) {
        fprintf( fp, "\"%s\" ", list.getName(i).c_str() );
    }
    fprintf( fp, "\n" );
    fprintf( fp, "ZONE T=\"Rectangular zone\"\n" );
    fprintf( fp, "I=%d, J=%d, K=1, ZONETYPE=Ordered\n", nx-1, ny-1 );
    fprintf( fp, "DATAPACKING=POINT\n" );
    fprintf( fp, "DT=(");
    for (int i=0; i<numVars; ++i) {
        fprintf( fp, "SINGLE ");
    }
    fprintf( fp, ")\n" );
    
    // Write the data
    const int lev=0; // finest grid
    for (int j=1; j<ny; ++j) {
        for (int i=1; i<nx; ++i) {
            for (int ind=0; ind < numVars; ++ind ) {
                fprintf( fp, "%.5e ", (*list.getVariable(ind))(lev,i,j) );
            }
            fprintf( fp, "\n" );
        }
    }

    fclose(fp);
    return true;
}

} // namespace ibpm
