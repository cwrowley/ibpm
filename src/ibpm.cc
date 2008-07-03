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
#include "Grid.h"
#include "Geometry.h"
#include "NonlinearNavierStokes.h"
#include "RungeKutta3.h"

using namespace std;

/*! \brief Sample main routine for IBFS code
 *  Set up a timestepper and advance the flow in time.
 *  Just a skeleton of this for now.  Eventually this should read parameters
 *  from an input file, do checking for bad input, etc.
 */
int main(int argc, char* argv[]) {
	cout << "Hello world!\n";
	
	// Setup grid
	int nx = 100;
	int ny = 100;
	double length = 2.0;
	Grid grid(nx, ny, length);
	
	// Setup geometry
	ifstream infile("geom.inp");
	assert(infile.good());
	Geometry geom;
	geom.load(infile);

	// Setup equations to solve
	double Reynolds=100;
	NonlinearNavierStokes model(grid, geom, Reynolds);
	model.init();

	// Setup timestepper
	dt = 0.01;
	RungeKutta3 solver(model, dt);

	// Load initial condition
	icFile = "initial.bin";
	State x;
	x.loadRestartFile(icFile);	

	// Step
	int numSteps = 100;
	for(int i=1; i <= numSteps; ++i) {
		cout << "step " << i << endl;
		solver.advance(x);		

		// Add something in here to output force data, probe data,
		// save restart files, etc.
		// Define a Logger class to handle all of this in one place?
	}
		
	return 0;
}