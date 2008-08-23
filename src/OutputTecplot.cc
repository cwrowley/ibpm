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
// $Revision:  $
// $LastChangedDate:  $
// $LastChangedBy: clancy $
// $HeadURL: $

#include "OutputTecplot.h"
#include "State.h"
#include "Output.h"
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;

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
    Grid grid = state.gamma.getGrid();
    int nx = grid.getNx();
    int ny = grid.getNy();
    double dx = grid.getDx();

    // Calculate the variables for output
    // Calculate the grid
    Scalar x(grid);
    Scalar y(grid);
    for (int i=0; i<=nx; ++i) {
        for (int j=0; j<=ny; ++j) {
            x(i,j) = grid.getXEdge(i);
            y(i,j) = grid.getYEdge(j);
        }
    }
    
    // Calculate vorticity
    Scalar vort = state.gamma;
    vort /= dx * dx;

    // Calculate velocities
    Scalar u(grid);
    Scalar v(grid);
    // TODO: Calculate u and v from fluxes
    for (int i=1; i<nx; ++i) {
        for (int j=1; j<ny; ++j) {
            u(i,j) = 0.5 * ( state.q(X,i,j) + state.q(X,i,j-1) ) / dx;
            v(i,j) = 0.5 * ( state.q(Y,i,j) + state.q(Y,i-1,j) ) / dx;
        }
    }
        
    // Store pointers to variables and corresponding names in vectors
    VarList list;
    list.addVariable( &x, "x" );
    list.addVariable( &y, "y" );
    list.addVariable( &u, "u" );
    list.addVariable( &v, "v" );
    list.addVariable( &vort, "Vorticity" );

    // Write the Tecplot file
    bool status = writeTecplotFileASCII( filename, title, list );
    return status;
}

// Write a Tecplot file with the specified filename and title,
// with data given in the VarList passed in
bool writeTecplotFileASCII(
    const char* filename,
    const char* title,
    const VarList& list
    ) {

    int numVars = list.getNumVars();
    assert( numVars > 0 );
    // Get grid information
    Grid grid = list.getVariable(0)->getGrid();
    int nx = grid.getNx();
    int ny = grid.getNy();

    // Write the header for the Tecplot file
    cout << "Writing Tecplot file " << filename << endl;
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
    for (int j=1; j<ny; ++j) {
        for (int i=1; i<nx; ++i) {
            for (int ind=0; ind < numVars; ++ind ) {
                fprintf( fp, "%.5e ", (*list.getVariable(ind))(i,j) );
            }
            fprintf( fp, "\n" );
        }
    }

    fclose(fp);
    return true;
}
