// Oseen.cc -- test Oseen vortex exact solution with IBPM code
//
// The Oseen vortex is an exact solution of Navier-Stokes, and has the form
//    \omega = \Gamma * Re / (4 \pi t) * exp( -r^2 Re / (4 t) )
//    u_\theta = \Gamma / (2 \pi r) * ( 1 - exp( -r^2 Re / (4 t)) )
//
// Here, the equations have been nondimensionalized by a length R and
// velocity U.  After Colonius and Taira (2007), we choose U to be the
// maximum velocity at the initial time t0, and choose R to be the radius
// at which this max velocity is attained.
//
// Clancy Rowley
// 22 Aug 2008
//
// $Revision$
// $Date$
// $Author$
// $HeadURL$

#include <iostream>
#include <sys/stat.h>
#include <math.h>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

void initializeOseenVortex(
    double Reynolds,
    State& x
);

void computeExactSolution(
    double Reynolds,
    double time,
    int timestep,
    State& exact
);

const double PI = 4. * atan(1.);

// Alpha is the solution to
//    e^{-\alpha} = 1 / (1 + 2\alpha)
// Setting
//    \omega = \Gamma * \alpha / PI * exp( -\alpha r^2 )
// makes the max velocity occur at r = 1
// and so corresponds by nondimensionalizing by the radius of u_max
const double alpha = 1.256431208626170;


int main(int argc, char* argv[]) {
	cout << "Test of Oseen vortex\n";
	
	// Setup grid
	int nx = 100;
	int ny = 100;
    int ngrid = 1;
	double length = 10;
    double xOffset = -5;
    double yOffset = -5;
	Grid grid( nx, ny, ngrid, length, xOffset, yOffset );
	
	// Empty geometry -- no body
	Geometry geom;
    int numPoints = 0;

	// Setup equations to solve
	double Reynolds=100;
    // No background flow
    BaseFlow q0( grid );

	NavierStokesModel model( grid, geom, Reynolds, q0 );
    model.init();
    
	// Setup timestepper
	double dt = 0.05;
    NonlinearIBSolver solver( grid, model, dt, Scheme::EULER );
    solver.init();

	// Setup initial condition
	State x( grid, numPoints );
    initializeOseenVortex( Reynolds, x );
    State exact( grid, numPoints );
    initializeOseenVortex( Reynolds, exact );
    State error( grid, numPoints );

    // Create output directory, if does not already exist
    mkdir( "oseen_out", S_IRWXU | S_IRWXG | S_IRWXO );
    
    // Setup output routines
    OutputTecplot outputComputed( "oseen_out/ibpm%03d.plt",
        "Oseen vortex, numerical, step %03d", 0 );
    OutputTecplot outputExact( "oseen_out/exact%03d.plt",
        "Oseen vortex, exact, step %03d", 0 );
    OutputTecplot outputError( "oseen_out/error%03d.plt",
        "Oseen vortex, error, step %03d", 0 );
    
    // Output initial condition
    outputComputed.doOutput( x );
    outputExact.doOutput( exact );
    outputError.doOutput( error );

	// Step
	int numSteps = 200;
    int iskip = 20;
	for(int i=1; i <= numSteps; ++i) {
		cout << "step " << i << endl;
		solver.advance( x );

        computeExactSolution( Reynolds, x.time, x.timestep, exact );

        // Compute error
        error.omega = x.omega - exact.omega;
        error.q = x.q - exact.q;
        error.time = x.time;
        error.timestep = x.timestep;

        // Write to output files
        if (i % iskip == 0) {
            outputComputed.doOutput( x );
            outputExact.doOutput( exact );
            outputError.doOutput( error );
        }
	}
 	return 0;
}

void computeExactSolution(
    double Reynolds,
    double time,
    int timestep,
    State& exact
    ) {
    const Grid& grid = exact.omega.getGrid();
    int nx = grid.Nx();
    int ny = grid.Ny();
    double dx = grid.Dx();

    const double circulation = 2 * PI * (1 + 1./(2*alpha) );
    double amp = circulation * Reynolds / (4 * PI * time);

    exact.time = time;
    exact.timestep = timestep;

    // Vorticity
    for (int i=1; i < nx; ++i ){
        double x = grid.getXEdge(0,i);
        for (int j=1; j < ny; ++j ){
            double y = grid.getYEdge(0,j);
            double rSquared = x*x + y*y;
            exact.omega(0,i,j) = amp * exp( -rSquared * Reynolds / (4 * time) );
        }
    }

    // Flux
    Flux::index ind;
    // X-direction flux
    for (ind = exact.q.begin(X); ind != exact.q.end(X); ++ind) {
        double x = exact.q.x(0,ind);
        double y = exact.q.y(0,ind);
        double r = sqrt( x*x + y*y );
        double uTheta = circulation / (2 * PI * r);
        uTheta *= 1 - exp( -r*r * Reynolds / (4 * time) );
        double theta = atan2( y, x );
        double u = -uTheta * sin(theta);
        exact.q(0,ind) = u * dx;
    }
    // Y-direction flux
    for (ind = exact.q.begin(Y); ind != exact.q.end(Y); ++ind) {
        double x = exact.q.x(0,ind);
        double y = exact.q.y(0,ind);
        double r = sqrt( x*x + y*y );
        double uTheta = circulation / (2 * PI * r);
        uTheta *= 1 - exp( -r*r * Reynolds / (4 * time) );
        double theta = atan2( y, x );
        double v =  uTheta * cos(theta);
        exact.q(0,ind) = v * dx;
    }
}

void initializeOseenVortex(
    double Reynolds,
    State& x
    ) {
    
    double t0 = Reynolds / (4 * alpha);
    computeExactSolution( Reynolds, t0, 0, x );    
}
