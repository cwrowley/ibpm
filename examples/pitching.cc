// plunging - example of main routine for a pitching flat plate
//
// Clancy Rowley
// 30 Aug 2008
// Princeton University
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

int main(int argc, char* argv[]) {
    cout << "Pitching flat plate example\n";

    // lift and drag
    double lift = 0.;
    double drag = 0.;

    // Setup grid
    int nx = 100;
    int ny = 100;
    int ngrid = 1;
    double length = 4;
    double xOffset = -1;
    double yOffset = -2;
    Grid grid( nx, ny, ngrid, length, xOffset, yOffset );

    // Make a flat plate, length 1, with center at 1/4 chord
    RigidBody plate;
    plate.addLine( 0, 0, 1, 0, grid.Dx() );
    plate.setCenter( 0.25, 0 );

    // Set the motion to pitching, amplitude = 0.25, period 10 time units
    double amplitude = 0.25;
    double freq = 0.1;
    PitchPlunge motion( amplitude, freq, 0, 0 );
    plate.setMotion( motion );
    Geometry geom;
    geom.addBody( plate );
    geom.moveBodies( 0 );

    // Setup equations to solve
    double Reynolds = 100;
    double magnitude = 1;
    double alpha = 0;  // angle of background flow
    BaseFlow q_potential( grid, magnitude, alpha );
    cout << "Setting up Navier Stokes model..." << flush;
    NavierStokesModel model( grid, geom, Reynolds, q_potential );
    model.init();
    cout << "done" << endl;

    // Setup timestepper
    double dt = 0.005;
    NonlinearIBSolver solver( grid, model, dt, Scheme::AB2 );
    solver.init();

    // Build the state variable, zero initial conditions
    State x(grid, geom.getNumPoints());
    x.omega = 0.;
    x.f = 0.;
    x.q = 0.;

    // Create output directory, if does not already exist
    mkdir( "pitching_out", S_IRWXU | S_IRWXG | S_IRWXO );

    // Setup output routines
    OutputForce force( "pitching_out/force.dat" );
    OutputTecplot tecplot( "pitching_out/pitch%03d.plt", "Pitching plate, step %03d" );
    Logger logger;
    // Output Tecplot file every few timesteps
    logger.addOutput( &tecplot, 25 );
    logger.addOutput( &force, 1 );
    logger.init();
    logger.doOutput( x );

    // Step
    const double PI = 4. * atan(1.);
    int numSteps = 250;
    for(int i=1; i <= numSteps; ++i) {
        double theta = amplitude * sin( 2 * PI * freq * x.time );
        cout << "step " << setw(4) << i
            << "  time = " << setw(5) << x.time
            << "  theta = " << theta << endl;
        solver.advance( x );
        x.computeNetForce( drag, lift );
        cout << " x force : " << setw(16) << drag*2 << " , y force : "
            << setw(16) << lift*2 << "\n";
        logger.doOutput( x );
    }
    logger.cleanup();
    return 0;
}
