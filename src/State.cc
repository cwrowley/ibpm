// State.cc
//
// Description:
// Implementation of the State class
//
// Author(s):
// Clancy Rowley
//
// Date: 1 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "State.h"
#include <string>
#include <stdio.h>

using namespace ibpm;

namespace ibpm {

State::State() {
    timestep = 0;
    time = 0.;
}

State::State(const Grid& grid, int numPoints ) {
    timestep = 0;
    time = 0.;
    resize( grid, numPoints );
}

void State::resize( const Grid& grid, int numPoints ) {
    q.resize( grid );
    gamma.resize( grid );
    f.resize( numPoints );
}

State::State( string filename ) {
    timestep = 0;
    time = 0.;
    load( filename );
}

State::~State() {}

bool State::load(const std::string& filename) {

    cerr << "Reading restart file " << filename << "..." << flush;
    FILE* fp = fopen( filename.c_str(), "rb" );
    if ( fp == NULL ) return false;

    // read Grid info
    int nx;
    int ny;
    double dx;
    double x0;
    double y0;
    
    fread( &nx, sizeof( int ), 1, fp );
    fread( &ny, sizeof( int ), 1, fp );
    fread( &dx, sizeof( double ), 1, fp );
    fread( &x0, sizeof( double ), 1, fp );
    fread( &y0, sizeof( double ), 1, fp );

    int numPoints;
    // read Geometry info
    fread( &numPoints, sizeof( int ), 1, fp );
    
    // check that Grid and Geometry in file match those expected
    bool success = true;
    if ( nx != q.Nx() || 
        ny != q.Ny() ||
        dx != q.Dx() ||
        x0 != q.getXEdge(0) ||
        y0 != q.getYEdge(0) ||
        numPoints != f.getNumPoints() ) {
        
        // If old grid was previously allocated, print a warning and set
        // the return value to false
        if ( q.Nx() > 0 ) {
            cerr << "Warning: grids do not match.  Resizing grid." << endl;
            success = false;
        }
        Grid newgrid( nx, ny, 1, dx * nx, x0, y0 );
        resize( newgrid, numPoints );
    }

    // read Flux q
    Flux::index qind;
    for ( qind = q.begin(); qind != q.end(); ++qind ) {
        success = success && fread( &( q(qind) ), sizeof( double ), 1, fp );
    }

    // read Scalar gamma
    for (int i=0; i<=nx; ++i ) {
        for ( int j=0; j<=ny; ++j ) {
            success = success && fread( &( gamma(i,j) ), sizeof( double ), 1, fp );
        }
    }

    // read BoundaryVector f
    for ( int i=0; i < numPoints; ++i ) {
        success = success && fread( &( f(X,i) ), sizeof( double ), 1, fp );
        success = success && fread( &( f(Y,i) ), sizeof( double ), 1, fp );        
    }

    // read timestep and time
    success = success && fread( &timestep, sizeof(int), 1, fp );
    success = success && fread( &time, sizeof(double), 1, fp );

    // close file
    fclose( fp );
    cerr << "done" << endl;
    return success;
}

bool State::save(std::string filename) const {
    cerr << "Writing restart file " << filename << "..." << flush;
    // open file
    FILE* fp = fopen( filename.c_str(), "wb" );
    if ( fp == NULL ) return false;

    // write Grid info
    const Grid& grid = q.getGrid();
    int nx = grid.Nx();
    int ny = grid.Ny();
    double dx = grid.Dx();
    double x0 = grid.getXEdge(0);
    double y0 = grid.getYEdge(0);
    
    fwrite( &nx, sizeof( int ), 1, fp );
    fwrite( &ny, sizeof( int ), 1, fp );
    fwrite( &dx, sizeof( double ), 1, fp );
    fwrite( &x0, sizeof( double ), 1, fp );
    fwrite( &y0, sizeof( double ), 1, fp );
        
    // write Geometry info
    int numPoints = f.getNumPoints();
    fwrite( &numPoints, sizeof( int ), 1, fp );
        
    // write Flux q
    Flux::index qind;
    for ( qind = q.begin(); qind != q.end(); ++qind ) {
        double qq = q(qind);
        fwrite( &qq, sizeof( double ), 1, fp );
    }

    // write Scalar gamma
    for (int i=0; i<=nx; ++i ) {
        for ( int j=0; j<=ny; ++j ) {
            double g = gamma(i,j);
            fwrite( &g, sizeof( double ), 1, fp );
        }
    }

    // write BoundaryVector f
    for ( int i=0; i < numPoints; ++i ) {
        double fx = f(X,i);
        double fy = f(Y,i);
        fwrite( &fx, sizeof( double ), 1, fp );
        fwrite( &fy, sizeof( double ), 1, fp );        
    }

    // write timestep and time
    fwrite( &timestep, sizeof(int), 1, fp );
    fwrite( &time, sizeof(double), 1, fp );

    // close file
    fclose( fp );
    cerr << "done" << endl;
    return true;
}

} // namespace ibpm
