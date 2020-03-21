// checkgeom - read in geometry file and write corresponding Tecplot file
//
// Use this utility to check the spacing of the boundary points.
// If the boundary points are too sparse, the boundary will look "leaky".
//
// NOTE: If the points are too finely spaced, the solution will not converge
// when used in a timestepper.
//
// Clancy Rowley
// 18 May 2009
// Princeton University
//
// $Revision:  $
// $LastChangedDate:  $
// $LastChangedBy: cwrowley $
// $HeadURL:  $

#include <iostream>
#include <iomanip>
#include "ibpm.h"
#include "Regularizer.h"

using namespace ibpm;

int main(int argc, char* argv[]) {
    cout << "Check geometry\n";

    ParmParser parser( argc, argv );
    bool helpFlag = parser.getFlag( "h", "print this help message and exit" );
    int nx = parser.getInt(
        "nx", "number of gridpoints in x-direction", 200 );
    int ny = parser.getInt(
        "ny", "number of gridpoints in y-direction", 200 );
    int ngrid = parser.getInt(
        "ngrid", "number of grid levels for multi-domain scheme", 1 );
    double length = parser.getDouble(
        "length", "length of finest domain in x-dir", 4.0 );
    double xOffset = parser.getDouble(
        "xoffset", "x-coordinate of left edge of finest domain", -2. );
    double yOffset = parser.getDouble(
        "yoffset", "y-coordinate of bottom edge of finest domain", -2. );
    string geomFile = parser.getString(
        "geom", "filename for reading geometry", "ibpm.geom" );
    string outFileName = parser.getString(
        "o", "filename for writing Tecplot file", "" );
    
    if ( ! parser.inputIsValid() || helpFlag ) {
        parser.printUsage( cerr );
        exit(1);
    }
    
    // Setup grid
    cout << "Grid parameters:" << endl
        << "  nx      " << nx << endl
        << "  ny      " << ny << endl
        << "  ngrid   " << ngrid << endl
        << "  length  " << length << endl
        << "  xoffset " << xOffset << endl
        << "  yoffset " << yOffset << endl;
    Grid grid( nx, ny, ngrid, length, xOffset, yOffset );

    // Setup geometry
    Geometry geom;
    cout << "Reading geometry from file " << geomFile << endl;
    if ( geom.load( geomFile ) ) {
        cout << "  " << geom.getNumPoints() << " points on the boundary" << endl;
    }
    else {
        cout << "  There were errors reading the geometry file." << endl;
        return 1;
    }
    
    if ( outFileName == "" ) return 0;

    Regularizer regularizer( grid, geom );
    regularizer.update();
    State x( grid, geom.getNumPoints() );

    // set the boundary force to 1 (in x- and y- directions)
    x.f = 1.;
    x.q = regularizer.toFlux( x.f );
    
    OutputTecplot tecplot( outFileName, "Check geometry", false );
    BaseFlow q_potential(grid, 0., 0.);
    tecplot.doOutput(q_potential, x);

    return 0;
}
