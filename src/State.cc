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
    
State::State(const Grid& grid, const Geometry& geom) :
    q(grid),
    gamma(grid),
    f( geom.getNumPoints() ),
    timestep(0),
    time(0.)
{}

State::~State() {}

#define CHECK(a,b,c)                                            \
    if ( !( (a) == (b) ) ) {                                    \
        cerr << "Error: " << c << " doesn't match." << endl     \
            << "  expected: " << (b) << endl                    \
            << "     read: " << (a) << endl;                    \
        return false;                                           \
    }

bool State::load(const std::string& filename) {

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
    const Grid& grid = q.getGrid();
    CHECK( nx, grid.getNx(), "nx" );
    CHECK( ny, grid.getNy(), "ny" );
    CHECK( dx, grid.getDx(), "dx" );
    CHECK( x0, grid.getXEdge(0), "x0" );
    CHECK( y0, grid.getYEdge(0), "x0" );

    CHECK( numPoints, f.getNumPoints(), "numPoints" );

    // read Flux q
    Flux::index qind;
    for ( qind = q.begin(); qind != q.end(); ++qind ) {
        fread( &( q(qind) ), sizeof( double ), 1, fp );
    }

    // read Scalar gamma
    for (int i=0; i<=nx; ++i ) {
        for ( int j=0; j<=ny; ++j ) {
            fread( &( gamma(i,j) ), sizeof( double ), 1, fp );
        }
    }

    // read BoundaryVector f
    for ( int i=0; i < numPoints; ++i ) {
        fread( &( f(X,i) ), sizeof( double ), 1, fp );
        fread( &( f(Y,i) ), sizeof( double ), 1, fp );        
    }

    // read timestep and time
    fread( &timestep, sizeof(int), 1, fp );
    fread( &time, sizeof(double), 1, fp );

    // close file
    fclose( fp );
    return true;
}

bool State::save(std::string filename) const {
    cerr << "Writing restart file " << filename << "..." << flush;
    // open file
    FILE* fp = fopen( filename.c_str(), "wb" );
    if ( fp == NULL ) return false;

    // write Grid info
    const Grid& grid = q.getGrid();
    int nx = grid.getNx();
    int ny = grid.getNy();
    double dx = grid.getDx();
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
