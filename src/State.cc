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
    omega.resize( grid );
    f.resize( numPoints );
}

State::State( string filename ) {
    timestep = 0;
    time = 0.;
    load( filename );
}

State::~State() {}
    
// Routine for computing X & Y forces
// Note that f is actually a body force (force per unit mass), so the net
// force is the integral over the domain.  By a property of the discrete
// delta function, this equals a sum of the BoundaryVector values times dx^2
void State::computeNetForce( double& xforce, double& yforce) const {
    xforce = 0;
    yforce = 0;
    for( int i=0; i<f.getNumPoints(); i++ ) {
        xforce += f(X,i);
        yforce += f(Y,i);
    }
    double dx2 = omega.Dx() * omega.Dx();
    xforce *= dx2;
    yforce *= dx2;
}

bool State::load(const std::string& filename) {

    cerr << "Reading restart file " << filename << "..." << flush;
    FILE* fp = fopen( filename.c_str(), "rb" );
    if ( fp == NULL ) return false;

    // read Grid info
    int nx;
    int ny;
    int ngrid;
    double dx;
    double x0;
    double y0;
    //double xShift;
    //double yShift;
    
    fread( &nx, sizeof( int ), 1, fp );
    fread( &ny, sizeof( int ), 1, fp );
    fread( &ngrid, sizeof( int ), 1, fp );
    fread( &dx, sizeof( double ), 1, fp );
    fread( &x0, sizeof( double ), 1, fp );
    fread( &y0, sizeof( double ), 1, fp );
    //fread( &xShift, sizeof( double ), 1, fp );
    //fread( &yShift, sizeof( double ), 1, fp );
    
    int numPoints;
    // read Geometry info
    fread( &numPoints, sizeof( int ), 1, fp );

    // check that Grid and Geometry in file match those expected
    bool success = true;
    if ( nx != q.Nx() || 
        ny != q.Ny() ||
        ngrid != q.Ngrid() ||
        dx != q.Dx() ||
        x0 != q.getXEdge(0,0) ||
        y0 != q.getYEdge(0,0) ||
        //xShift != q.getXShift() ||
        //yShift != q.getYShift() ||
        numPoints != f.getNumPoints() ) {
        
        // If old grid was previously allocated, print a warning and set
        // the return value to false
        if ( q.Nx() > 0 ) {
            cerr << "Warning: grids do not match.  Resizing grid." << endl;
            /*cerr << "    nx:        " << ( nx == q.Nx() ) << endl
                 << "    ny:        " << ( ny == q.Ny() ) << endl
                 << "    ngrid:     " << ( ngrid == q.Ngrid() ) << endl
                 << "    dx:        " << ( dx == q.Dx() ) << endl
                 << "    x0:        " << ( x0 == q.getXEdge(0,0) ) << endl
                 << "    y0:        " << ( y0 == q.getYEdge(0,0) ) << endl
                 << "    numPoints: " << ( numPoints == f.getNumPoints() ) << endl
                 << endl;*/
            success = false;
        }
        //Grid newgrid( nx, ny, ngrid, dx * nx, x0, y0, xShift, yShift );
        Grid newgrid( nx, ny, ngrid, dx * nx, x0, y0 );
        resize( newgrid, numPoints );
    }

    // read Flux q
    Flux::index qind;
    for ( int lev=0; lev < q.Ngrid(); ++lev ) {
        for ( qind = q.begin(); qind != q.end(); ++qind ) {
            success = success && fread( &( q(lev,qind) ), sizeof( double ), 1, fp );
        }
    }
    
    // read Scalar omega
    for ( int lev=0; lev < q.Ngrid(); ++lev ) {
        for (int i=1; i<nx; ++i ) {
            for ( int j=1; j<ny; ++j ) {
                success = success && fread( &( omega(lev,i,j) ), sizeof( double ), 1, fp );
            }
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
    int ngrid = grid.Ngrid();
    double dx = grid.Dx();
    double x0 = grid.getXEdge(0,0);
    double y0 = grid.getYEdge(0,0);
    //double xShift = grid.getXShift();
    //double yShift = grid.getYShift();
    
    fwrite( &nx, sizeof( int ), 1, fp );
    fwrite( &ny, sizeof( int ), 1, fp );
    fwrite( &ngrid, sizeof( int ), 1, fp );
    fwrite( &dx, sizeof( double ), 1, fp );
    fwrite( &x0, sizeof( double ), 1, fp );
    fwrite( &y0, sizeof( double ), 1, fp );
    //fwrite( &xShift, sizeof( double ), 1, fp );
    //fwrite( &yShift, sizeof( double ), 1, fp );
        
    // write Geometry info
    int numPoints = f.getNumPoints();
    fwrite( &numPoints, sizeof( int ), 1, fp );
        
    // write Flux q
    Flux::index qind;
    for ( int lev=0; lev < q.Ngrid(); ++lev ) {
        for ( qind = q.begin(); qind != q.end(); ++qind ) {
            double qq = q(lev,qind);
            fwrite( &qq, sizeof( double ), 1, fp );
        }
    }

    // write Scalar omega
    for ( int lev=0; lev < q.Ngrid(); ++lev ) {
        for (int i=1; i<nx; ++i ) {
            for ( int j=1; j<ny; ++j ) {
                double g = omega(lev,i,j);
                fwrite( &g, sizeof( double ), 1, fp );
            }
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
