#include <iostream>
#include <sys/stat.h>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

int main(int argc, char* argv[]) {
    // Setup grid
    int nx = 100;
    int ny = 50;
    int ngrid = 3;
    double length = 10;
    double xOffset = -2;
    double yOffset = -2.5;
    Grid grid( nx, ny, ngrid, length, xOffset, yOffset );
    
    // Make a cylinder, diameter 1, with center at (0,0)
    RigidBody cylinder;
    cylinder.addCircle( 0, 0, 0.5, grid.Dx() );
    Geometry geom;
    geom.addBody(cylinder);

    // Setup equations to solve, with uniform flow U = 1
    double magnitude = 1.;
    double angle = 0.;
    BaseFlow q_potential( grid, magnitude, angle );
    double Reynolds=100;
    NavierStokesModel model( grid, geom, Reynolds, q_potential );
    model.init();

    // Setup timestepper
    double dt = 0.05;
    NonlinearIBSolver solver( grid, model, dt, Scheme::EULER );
    solver.init();

    // State variable with zero initial conditions
    State x( grid, geom.getNumPoints() );
    x.omega = 0.;

    // Setup output routines
    OutputTecplot tecplot( "cyl%03d.plt", "Cylinder, step %03d", false );
    
    // Output initial condition
    tecplot.doOutput( x );

    // Step
    int numSteps = 10;
    for(int i=1; i <= numSteps; ++i) {
        cout << "step " << i << endl;
        solver.advance( x );
    }

    // Output final state
    tecplot.doOutput( x );
    return 0;
}
