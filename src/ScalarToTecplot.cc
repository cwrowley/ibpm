
#include "ScalarToTecplot.h"
#include <stdio.h>
#include <cstring>
#include <vector>

namespace ibpm {
    
class VarList {
public:
    void addVariable( const Scalar* var, string varName ) {
        if ( getNumVars() > 0 ) {
           assert( _vars[0]->getGrid().isEqualTo( var->getGrid() ) );
        }
        _vars.push_back( var );
        _names.push_back( varName );
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
    
bool writeTecplotFileASCII( const char* filename, const char* title, const VarList& list, int lev ) {
    int numVars = list.getNumVars();
    assert( numVars > 0 );
    
    // Get grid information
    const Grid& grid = list.getVariable(0)->getGrid();
    assert( lev < grid.Ngrid() );
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
    
    
bool ScalarToTecplot( vector<const Scalar*> varVec, vector<string> varNameVec, string filename, string title, int lev ) {
    assert( varVec.size() > 0 );
    assert( varVec.size() == varNameVec.size() );
    
    // Get grid dimensions
    const Grid& grid = varVec[0]->getGrid();
    int nx = grid.Nx();
    int ny = grid.Ny();
    int ngrid = grid.Ngrid();
    assert( lev < ngrid );
    
    // Calculate the variables for output
    // Calculate the grid
    Scalar x(grid);
    Scalar y(grid);
    for (int _lev=0; _lev<ngrid; ++_lev) {
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny; ++j) {
                x(_lev,i,j) = grid.getXEdge(_lev,i);
                y(_lev,i,j) = grid.getYEdge(_lev,j);
            }
        }
    }
    
    // Store pointers to variables and corresponding names in vectors
    VarList list;
    list.addVariable( &x, "x" );
    list.addVariable( &y, "y" );
    for( unsigned int i = 0; i < varVec.size(); i++ ) {
        list.addVariable( varVec[i], varNameVec[i] );
    }
    
    // Add timestep to filename and title
    char _filename[BUFSIZ];
    strncpy( _filename, filename.c_str(), BUFSIZ-1 );
    char _title[BUFSIZ];
    strncpy( _title, title.c_str(), BUFSIZ-1 );
    
    // Write the Tecplot file
    bool status = writeTecplotFileASCII( _filename, _title, list, lev );
    return status;
}
    
bool ScalarToTecplot( const Scalar* var, string varName, string filename, string title, int lev ) {
    vector<const Scalar*> varVec;
    varVec.push_back( var );
    
    vector<string> varNameVec;
    varNameVec.push_back( varName );
    
    bool status = ScalarToTecplot( varVec, varNameVec, filename, title, lev);
    return status;
}
    
bool ScalarToTecplot( vector<const Scalar*> varVec, vector<string> varNameVec, string filename, string title ) { 
    return ScalarToTecplot( varVec, varNameVec, filename, title, 0 );
}
    
bool ScalarToTecplot( const Scalar* var, string varName, string filename, string title ) {
    return ScalarToTecplot( var, varName, filename, title, 0 );
}
} // namespace ibpm
