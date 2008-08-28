/*!
\file ibpm.cc

\brief Sample main routine for IBFS code

\author Clancy Rowley
\date  3 Jul 2008

$Revision$
$Date$
$Author$
$HeadURL$
*/

#include <iostream>
#include <fstream>
#include <string>
#include "Grid.h"
#include "RigidBody.h"
#include "Geometry.h"
#include "NavierStokesModel.h"
#include "Euler.h"
#include "State.h"
#include "Logger.h"
#include "OutputTecplot.h"
#include "VectorOperations.h"

using namespace std;

/*! \brief Sample main routine for IBFS code
 *  Set up a timestepper and advance the flow in time.
 *  Just a skeleton of this for now.  Eventually this should read parameters
 *  from an input file, do checking for bad input, etc.
 */
int main(int argc, char* argv[]) {
	cout << "Hello world!\n";

     // lift and drag
     double lift = 0.;
     double drag = 0.;
	
	// Setup grid
	int nx = 400;
	int ny = 400;
	double length = 4.0;
    double xOffset = -2;
    double yOffset = -2;
	Grid grid( nx, ny, length, xOffset, yOffset );
	
	// Setup geometry
    // ifstream infile("geom.inp");
    // assert(infile.good());
	Geometry geom;
    RigidBody circle;
    double radius = 0.5;
    double dTheta = grid.getDx() / radius;
    double pi = 4. * atan(1.);
    int numPoints = 2 * pi / dTheta + 1;
    circle.addCircle( 0, 0, radius, numPoints );
    geom.addBody( circle );
    // geom.load(infile);

	// Setup equations to solve
	double Reynolds=100;
    double magnitude = 1;
    double alpha = 0;  // angle of background flow
    Flux q_potential = Flux::UniformFlow( grid, magnitude, alpha );
	NonlinearNavierStokes model( grid, geom, Reynolds, q_potential );
    // model.init();

	// Setup timestepper
	double dt = 0.005;
	Euler solver(model, dt);

	// Load initial condition
	string icFile = "initial.bin";
	State x(grid, geom);
	x.loadRestartFile(icFile);	
    x.gamma = 0.;

    // Setup output routines
    OutputTecplot tecplot( "ibpm%03d.plt", "Test run, step %03d" );
    Logger logger;
    // Output Tecplot file every timestep
    logger.addOutput( &tecplot, 25 );
    logger.init();
    logger.doOutput( x );
    
	// Step
	int numSteps = 5000;
	for(int i=1; i <= numSteps; ++i) {
		cout << "step " << i << endl;
		solver.advance( x );
                computeNetForce(x.f, drag, lift);
                cout << "x force : " << setw(16) << drag/dt << " , y force : " << setw(16) << lift/dt << "\n";
        logger.doOutput( x );
	}
    logger.cleanup();
	return 0;
}
