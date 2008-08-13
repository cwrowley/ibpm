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
#include "Geometry.h"
#include "NavierStokesModel.h"
#include "Euler.h"
#include "State.h"

using namespace std;

/*! \brief Sample main routine for IBFS code
 *  Set up a timestepper and advance the flow in time.
 *  Just a skeleton of this for now.  Eventually this should read parameters
 *  from an input file, do checking for bad input, etc.
 */
int main(int argc, char* argv[]) {
	cout << "Hello world!\n";
	
	// Setup grid
	int nx = 20;
	int ny = 20;
	double length = 2.0;
    double xOffset = -1;
    double yOffset = -1;
	Grid grid( nx, ny, length, xOffset, yOffset );
	
	// Setup geometry
    ifstream infile("geom.inp");
    assert(infile.good());
	Geometry geom;
	geom.load(infile);

	// Setup equations to solve
	double Reynolds=100;
    double magnitude = 1;
    double alpha = 0;  // angle of background flow
    Flux q_potential = Flux::UniformFlow( grid, magnitude, alpha );
	NonlinearNavierStokes model( grid, geom, Reynolds, q_potential );
    // model.init();

	// Setup timestepper
	double dt = 0.01;
	Euler solver(model, dt);

	// Load initial condition
	string icFile = "initial.bin";
	State x(grid, geom);
	x.loadRestartFile(icFile);	

	// Step
	int numSteps = 10;
	for(int i=1; i <= numSteps; ++i) {
		cout << "step " << i << endl;
		solver.advance(x);
		// Add something in here to output force data, probe data,
		// save restart files, etc.
		// Define a Logger class to handle all of this in one place?
	}
		
	return 0;
}