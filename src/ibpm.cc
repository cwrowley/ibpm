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
#include <iomanip>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

enum ModelType { LINEAR, NONLINEAR, ADJOINT, INVALID };

// Return a solver of the appropriate type (e.g. Euler, RK2 )
TimeStepper* GetSolver(
    Grid& grid,
    NavierStokesModel& model,
    double dt,
    string solverType
);

// Return the type of model specified in the string modelName
ModelType str2model( string modelName );

/*! \brief Main routine for IBFS code
 *  Set up a timestepper and advance the flow in time.
 */
int main(int argc, char* argv[]) {
    cout << "Immersed Boundary Projection Method (IBPM), version "
        << IBPM_VERSION << endl;

    // Get parameters
    ParmParser parser( argc, argv );
    bool helpFlag = parser.getFlag( "h", "print this help message and exit" );
    string name = parser.getString( "name", "run name", "ibpm" );
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
        "geom", "filename for reading geometry", name + ".geom" );
    double Reynolds = parser.getDouble("Re", "Reynolds number", 100.);
    double dt = parser.getDouble( "dt", "timestep", 0.01 );
    string modelName = parser.getString(
        "model", "type of model (linear, nonlinear, adjoint)", "nonlinear" );
    string baseFlow = parser.getString(
        "baseflow", "base flow for linear/adjoint model", "" );
    string integratorType = parser.getString(
        "scheme", "timestepping scheme (euler,ab2,rk2,rk3)", "rk2" );
    string icFile = parser.getString( "ic", "initial condition filename", "");
    string outdir = parser.getString(
        "outdir", "directory for saving output", "." );
    int iTecplot = parser.getInt(
        "tecplot", "if >0, write a Tecplot file every n timesteps", 100);
    int iRestart = parser.getInt(
        "restart", "if >0, write a restart file every n timesteps", 100);
    int iForce = parser.getInt(
        "force", "if >0, write forces every n timesteps", 1);
    int numSteps = parser.getInt(
        "nsteps", "number of timesteps to compute", 250 );

    ModelType modelType = str2model( modelName );
    
    if ( ! parser.inputIsValid() || modelType == INVALID || helpFlag ) {
        parser.printUsage( cerr );
        exit(1);
    }
    
    if ( modelType != NONLINEAR && baseFlow == "" ) {
        cout << "ERROR: for linear or adjoint models, "
            "must specify a base flow" << endl;
        exit(1);
    }

    // create output directory if not already present
    AddSlashToPath( outdir );
    mkdir( outdir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO );

    // output command line arguments
    string cmd = parser.getParameters();
    cout << "Command:" << endl << cmd << endl;
    parser.saveParameters( outdir + name + ".cmd" );

    // Name of this run
    cout << "Run name: " << name << endl;

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
        exit(-1);
    }

    // Setup equations to solve
    cout << "Reynolds number = " << Reynolds << endl;
    double magnitude = 1;
    double alpha = 0;  // angle of background flow
    Flux q_potential = Flux::UniformFlow( grid, magnitude, alpha );
    cout << "Setting up Navier Stokes model..." << flush;
    NavierStokesModel* model = NULL;
    if ( modelType == NONLINEAR ) {
        model = new NonlinearNavierStokes( grid, geom, Reynolds, q_potential);
    }
    else {
        State x0( grid, geom.getNumPoints() );
        x0.load( baseFlow );
        if ( modelType == LINEAR ) {
            model = new LinearizedNavierStokes( grid, geom, Reynolds, x0 );
        }
        else if ( modelType == ADJOINT ){
            model = new AdjointNavierStokes( grid, geom, Reynolds, x0 );
        }
    }
    assert( model != NULL );
    model->init();
    cout << "done" << endl;

    // Setup timestepper
    TimeStepper* solver = GetSolver( grid, *model, dt, integratorType );
    cout << "Using " << solver->getName() << " timestepper" << endl;
    cout << "  dt = " << dt << endl;
    if ( ! solver->load( outdir + name ) ) {
        solver->init();
        solver->save( outdir + name );
    }

    // Load initial condition
    State x( grid, geom.getNumPoints() );
    x.gamma = 0.;
    x.f = 0.;
    x.q = 0.;
    if ( icFile != "" ) {
        cout << "Loading initial condition from file: " << icFile << endl;
        if ( ! x.load(icFile) ) {
            cout << "  (failed: using zero initial condition)" << endl;
        }
    }
    else {
        cout << "Using zero initial condition" << endl;
    }

    // Setup output routines
    OutputTecplot tecplot( outdir + name + "%03d.plt", "Test run, step %03d");
    OutputRestart restart( outdir + name + "%03d.bin" );
    OutputForce force( outdir + name + ".force" ); 

    Logger logger;
    // Output Tecplot file every timestep
    if ( iTecplot > 0 ) {
        cout << "Writing Tecplot file every " << iTecplot << " steps" << endl;
        logger.addOutput( &tecplot, iTecplot );
    }
    if ( iRestart > 0 ) {
        cout << "Writing restart file every " << iRestart << " steps" << endl;
        logger.addOutput( &restart, iRestart );
    }
    if ( iForce > 0 ) {
        cout << "Writing forces every " << iForce << " steps" << endl;
        logger.addOutput( &force, iForce );
    }
    logger.init();
    logger.doOutput( x );
    
    cout << "Integrating for " << numSteps << " steps" << endl;

    for(int i=1; i <= numSteps; ++i) {
        cout << "step " << i << endl;
        solver->advance( x );
        double lift;
        double drag;
        computeNetForce( x.f, drag, lift );
        cout << "x force : " << setw(16) << drag*2 << " , y force : "
            << setw(16) << lift*2 << "\n";
        logger.doOutput( x );
    }
    logger.cleanup();
    delete solver;
    delete model;
    return 0;
}

ModelType str2model( string modelName ) {
    ModelType type;
    MakeLowercase( modelName );
    if ( modelName == "nonlinear" ) {
        type = NONLINEAR;
    }
    else if ( modelName == "linear" ) {
        type = LINEAR;
    }
    else if ( modelName == "adjoint" ) {
        type = ADJOINT;
    }
    else {
        cerr << "Unrecognized model: " << modelName << endl;
        type = INVALID;
    }
    return type;
}

TimeStepper* GetSolver(
    Grid& grid,
    NavierStokesModel& model,
    double dt,
    string solverType
    ) {
    MakeLowercase( solverType );
    if ( solverType == "euler" ) {
        return new Euler( grid, model, dt );
    }
    else if ( solverType == "ab2" ) {
        return new AdamsBashforth( grid, model, dt );
    }
    else if ( solverType == "rk2" ) {
        return new RungeKutta2( grid, model, dt );
    }
    else if ( solverType == "rk3" ) {
        return new RungeKutta3( grid, model, dt );
    }
    else {
        cerr << "ERROR: unrecognized solver: " << solverType << endl;
        exit(1);
    }
}

