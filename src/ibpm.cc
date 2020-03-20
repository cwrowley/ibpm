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

using namespace ibpm;

enum ModelType { LINEAR, NONLINEAR, ADJOINT, LINEARPERIODIC, SFD, INVALID };

// Return the type of model specified in the string modelName
ModelType str2model( string modelName );

// Return the integration scheme specified in the string integratorType
Scheme::SchemeType str2scheme( string integratorType );

/*! \brief Main routine for IBFS code
 *  Set up a timestepper and advance the flow in time.
 */
int main(int argc, char* argv[]) {
    cout << "Immersed Boundary Projection Method (IBPM), version "
        << IBPM_VERSION << "\n" << endl;

    // Get parameters
    ParmParser parser( argc, argv );
    bool helpFlag = parser.getFlag( "h", "print this help message and exit" );
    
    // Output parameters
    string name = parser.getString( "name", "run name", "ibpm" );
    string outdir = parser.getString( "outdir", "directory for saving output", "." );
    int iTecplot = parser.getInt( "tecplot", "if >0, write a Tecplot file every n timesteps", 100);
    bool TecplotAllGrids = parser.getBool( "tecplotallgrids", "Tecplot output for all grids, or not", false );
    int iRestart = parser.getInt( "restart", "if >0, write a restart file every n timesteps", 100);
    int iForce = parser.getInt( "force", "if >0, write forces every n timesteps", 1);
    int iEnergy = parser.getInt( "energy", "if >0, write energy every n timesteps", 0);
    string numDigitInFileName = parser.getString( "numdigfilename", "number of digits for time representation in filename", "%05d");
    
    // Grid parameters
    int nx = parser.getInt( "nx", "number of gridpoints in x-direction", 200 );
    int ny = parser.getInt( "ny", "number of gridpoints in y-direction", 200 );
    int ngrid = parser.getInt( "ngrid", "number of grid levels for multi-domain scheme", 1 );
    double length = parser.getDouble( "length", "length of finest domain in x-dir", 4.0 );
    double xOffset = parser.getDouble( "xoffset", "x-coordinate of left edge of finest domain", -2. );
    double yOffset = parser.getDouble( "yoffset", "y-coordinate of bottom edge of finest domain", -2. );
    double xShift = parser.getDouble( "xshift", "percentage offset between grid levels in x-direction", 0. );
    double yShift = parser.getDouble( "yshift", "percentage offset between grid levels in y-direction", 0. );
    double alpha = parser.getDouble( "alpha", "angle of attack of base flow", 0.);    

    // Simulation parameters
    string geomFile = parser.getString( "geom", "filename for reading geometry", name + ".geom" );
    bool ubf = parser.getBool( "ubf", "Use unsteady base flow, or not", false );
    double Reynolds = parser.getDouble("Re", "Reynolds number", 100.);
    string modelName = parser.getString( "model", "type of model (linear, nonlinear, adjoint, linearperiodic, sfd)", "nonlinear" );
    string baseFlow = parser.getString( "baseflow", "base flow for linear/adjoint model", "" );
    
    // Initial condition
    string icFile = parser.getString( "ic", "initial condition filename", "");
    bool resetTime = parser.getBool( "resettime", "Reset time when subtracting ic by baseflow (1/0(true/false))", false);
    bool subtractBaseflow = parser.getBool( "subbaseflow", "Subtract ic by baseflow (1/0(true/false))", false);
    
    // Integration parameters
    double dt = parser.getDouble( "dt", "timestep", 0.02 );
    int numSteps = parser.getInt( "nsteps", "number of timesteps to compute", 250 );
    string integratorType = parser.getString( "scheme", "timestepping scheme (euler,ab2,rk3,rk3b)", "rk3" );
    
    // Linear-periodic model
    int period = parser.getInt( "period", "period of periodic baseflow", 1);
    int periodStart = parser.getInt( "periodstart", "start time of periodic baseflow", 0);
    string periodBaseFlowName = parser.getString( "pbaseflowname", "name of periodic baseflow, e.g. 'flow/ibpmperiodic%05d.bin', with '%05d' as time, decided by periodstart/period", "" );
    
    // SFD
    double chi = parser.getDouble( "chi", "sfd gain", 0.02 );
    double Delta = parser.getDouble( "Delta", "sfd cutoff frequency", 15. );

    
    ModelType modelType = str2model( modelName );
    Scheme::SchemeType schemeType = str2scheme( integratorType );
    
    if ( ! parser.inputIsValid() || modelType == INVALID || helpFlag ) {
        parser.printUsage( cerr );
        exit(1);
    }
    
    // modify this long if statement?
    if ( ( modelType != NONLINEAR ) && ( modelType != SFD ) ) {
        if (modelType != LINEARPERIODIC && baseFlow == "" ){
            cout << "ERROR: for linear or adjoint models, "
            "must specify a base flow" << endl;
            exit(1);
        }
        else if (modelType != LINEARPERIODIC && periodBaseFlowName != ""){
            cout << "WARNING: for linear or adjoint models, "
            "a periodic base flow is not needed" << endl;
            exit(1);
        }
        else if (modelType == LINEARPERIODIC && periodBaseFlowName == "" ) {
            cout << "ERROR: for linear periodic model, "
            "must specify a periodic base flow" << endl;
            exit(1);
        }
        else if (modelType == LINEARPERIODIC && baseFlow != "" ) {
            cout << "WARNING: for linear periodic model, "
            "a single baseflow is not needed" << endl;
            exit(1);
        }       
    }
    
    // create output directory if not already present
    AddSlashToPath( outdir );
    mkdir( outdir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO );

    // output command line arguments
    string cmd = parser.getParameters();
    cout << "Command:" << endl << cmd << "\n" << endl;
    parser.saveParameters( outdir + name + ".cmd" );

    // Name of this run
    cout << "Run name: " << name << "\n" << endl;

    // Setup grid
    cout << "Grid parameters:" << endl
        << "    nx      " << nx << endl
        << "    ny      " << ny << endl
        << "    ngrid   " << ngrid << endl
        << "    length  " << length << endl
        << "    xoffset " << xOffset << endl
        << "    yoffset " << yOffset << endl
        << "    xshift  " << xShift << endl
        << "    yshift  " << yShift << endl
        << endl;
    Grid grid( nx, ny, ngrid, length, xOffset, yOffset, xShift, yShift );

    // Setup geometry
    Geometry geom;
    cout << "Reading geometry from file " << geomFile << endl;
    if ( geom.load( geomFile ) ) {
        cout << "    " << geom.getNumPoints() << " points on the boundary" << "\n" << endl;
    }
    else {
        exit(-1);
    }
    
    // Setup equations to solve
    cout << "Reynolds number = " << Reynolds << "\n" << endl;
    cout << "Setting up Immersed Boundary Solver..." << flush;
    double magnitude = 1;
//    double alpha = 0;  // angle of background flow
    double pi = 4. * atan(1.);
    alpha = alpha*pi/180.;
    double xC = 0, yC = 0;
    BaseFlow q_potential( grid, magnitude, alpha );
    // See if unsteady base flow can be used.  Only implemented for a single RigidBody in motion.  
    // In the future, have a function geom.ubfEligible() that will make sure that the first RigidBody is moving.
    if( ! geom.isStationary() && (geom.getNumBodies() == 1) && ubf ) {
        Motion* m = geom.transferMotion(); // pull motion from first RigidBody object
        geom.transferCenter(xC,yC);        // pull center of motion from RigidBody object
        q_potential.setMotion( *m ); 
        q_potential.setCenter(xC,yC);
    }
    if( ubf && (geom.getNumBodies() != 1) ) {
        cout << "Unsteady base flow is only supported for a single moving body.  Exiting program." << endl;
        exit(1);
    }
  
    NavierStokesModel* model = NULL;
    IBSolver* solver = NULL;
    SFDSolver* SFDsolver = NULL;
    State x00( grid, geom.getNumPoints() ); 

    switch (modelType){
        case NONLINEAR: 
            model =  new NavierStokesModel( grid, geom, Reynolds, q_potential );
            solver = new NonlinearIBSolver( grid, *model, dt, schemeType );
            break;
        case LINEAR:
            if ( ! x00.load( baseFlow ) ) {
                cout << "baseflow failed to load.  Exiting program." << endl;
                exit(1);
            }
            model =  new NavierStokesModel( grid, geom, Reynolds );
            solver = new LinearizedIBSolver( grid, *model, dt, schemeType, x00 );
            break;          
        case ADJOINT:
            if ( ! x00.load( baseFlow ) ) {
                cout << "baseflow failed to load.  Exiting program." << endl;
                exit(1);
            }
            model =  new NavierStokesModel( grid, geom, Reynolds );
            solver = new AdjointIBSolver( grid, *model, dt, schemeType, x00 );
            break;          
        case LINEARPERIODIC:{
            // load periodic baseflow files
            vector<State> x0(period, x00);  
            char pbffilename[256];
            string pbf = periodBaseFlowName;
            for (int i=0; i < period; i++) {
                //cout << "loading the " << i << "-th periodic baseflow:" << endl;
                sprintf(pbffilename, pbf.c_str(), i + periodStart);
                x0[i].load(pbffilename);                      
            }
            x00 = x0[0];
            model =  new NavierStokesModel( grid, geom, Reynolds );
            solver = new LinearizedPeriodicIBSolver( grid, *model, dt, schemeType, x0, period ) ;   
            break;
            }
        case SFD:{ 
            cout << "SFD parameters:" << endl;
            cout << "    chi =   " << chi << endl;
            cout << "    Delta = " << Delta << endl << endl;
            model =  new NavierStokesModel( grid, geom, Reynolds, q_potential );
            SFDsolver = new SFDSolver( grid, *model, dt, schemeType, Delta, chi ) ;
            solver = SFDsolver;
            break;
            }
        case INVALID:
            cout << "ERROR: must specify a valid modelType" << endl;
            exit(1);
            break;
    }
    
    assert( model != NULL );
    assert( solver != NULL );
    if( modelType == SFD ) {
        assert( chi != 0 );
        assert( SFDsolver != NULL );
    }
    // NOTE: still need to initialize model, but wait until after loading the initial
    //       condition, so we know what the initial time is, for moving the bodies

    // Load initial condition
    State x( grid, geom.getNumPoints() );
    x.omega = 0.;
    x.f = 0.;
    x.q = 0.;
    if (icFile != "") { 
        cout << "Loading initial condition from file: " << icFile << endl;
        if ( ! x.load(icFile) ) {
            cout << "    (failed: using zero initial condition)" << endl;
        }
        if ( subtractBaseflow == true ) {
            cout << "    Subtracting initial condition by baseflow to form a linear initial perturbation" << endl;
            if (modelType != NONLINEAR) {
                assert((x.q).Ngrid() == (x00.q).Ngrid());
                assert((x.omega).Ngrid() == (x00.omega).Ngrid());
                x.q -= x00.q;
                x.omega -= x00.omega;
                x.f = 0;
            }
            else {
                cout << "Flag subbaseflow should be true only for linear cases"<< endl;
                exit(1);
            }
        }
        
        if ( modelType == SFD ) {
            SFDsolver->loadFilteredState( icFile );
        }
         
    }
    else {
        cout << "Using zero initial condition" << endl;
    }
    
    if (resetTime == true) {
        x.timestep = 0;
        x.time = 0.;
    }
    
    // update the geometry to the current time
    geom.moveBodies( x.time );
    
    // Initialize model and timestepper
    model->init();
    cout << "using " << solver->getName() << " timestepper" << endl;
    cout << "    dt = " << dt << "\n" << endl;
    if ( ! solver->load( outdir + name ) ) {
        // Set the tolerance for a ConjugateGradient solver below
        // Otherwise default is tol = 1e-7
        // solver->setTol( 1e-8 )
        solver->init();
        solver->save( outdir + name );
    }
    
    // Calculate flux for state, in case only vorticity was saved
    if( ! q_potential.isStationary() ) {
       q_potential.setAlphaMag( x.time );
       alpha = q_potential.getAlpha();
    }
    model->updateOperators( x.time );
    model->refreshState( x );

    cout << endl << "Initial timestep = " << x.timestep << "\n" << endl;

    // Setup output routines
    OutputTecplot tecplot( outdir + name + numDigitInFileName + ".plt", "Test run, step" +  numDigitInFileName, TecplotAllGrids);
    if(TecplotAllGrids) tecplot.setFilename( outdir + name + numDigitInFileName + "_g%01d.plt" );
    OutputRestart restart( outdir + name + numDigitInFileName + ".bin" );
    OutputForce force( outdir + name + ".force" ); 
    OutputEnergy energy( outdir + name + ".energy" ); 
    
    Logger logger;
    // Output Tecplot file every timestep
    if ( iTecplot > 0 ) {
        cout << "Writing Tecplot file every " << iTecplot << " step(s)" << endl;
        logger.addOutput( &tecplot, iTecplot );
    }
    if ( iRestart > 0 ) {
        cout << "Writing restart file every " << iRestart << " step(s)" << endl;
        logger.addOutput( &restart, iRestart );
    }
    if ( iForce > 0 ) {
        cout << "Writing forces every " << iForce << " step(s)" << endl;
        logger.addOutput( &force, iForce );
    }
    if ( iEnergy > 0 ) {
        cout << "Writing energy every " << iForce << " step(s)" << endl;
        logger.addOutput( &energy, iEnergy );
    }
    cout << endl;
    logger.init();
    logger.doOutput( q_potential, x );

    cout << "Integrating for " << numSteps << " steps" << endl;
    for(int i=1; i <= numSteps; ++i) {
        cout << "\nstep " << i << endl; 
        State xtemp( x ); // For SFD norm calculation
        solver->advance( x );
        double lift;
        double drag;
        double xF, yF; // forces in x and y direction (same as drag,lift if alpha=0)
        x.computeNetForce( xF, yF );
        // If there is an unsteady base flow, transform body frame normal and parallel forces into lab frame lift and drag
        if( ! q_potential.isStationary() ) {
            q_potential.setAlphaMag(x.time);
            alpha = q_potential.getAlpha();
        }
        drag = xF * cos(alpha) + yF * sin(alpha);
        lift = xF * -1.*sin(alpha) + yF * cos(alpha);
        cout << "    x force: " << setw(16) << drag*2 << ", y force: "
            << setw(16) << lift*2 << "\n";
        logger.doOutput( q_potential, x );

        // For SFD
        if( modelType == SFD ) {
            // Inner product of fluxes is equal to inner product of vorticity (with weighted inner product for latter)
            Flux dq = xtemp.q-x.q;
            double q = sqrt( InnerProduct( x.q, x.q ) );
            double twoNorm = sqrt( InnerProduct( dq, dq ) ) / ( q * dt );
            
            if ( (x.timestep % iRestart == 0 ) && (chi != 0.0) ) {
                SFDsolver->saveFilteredState( outdir, name, numDigitInFileName );
            }
            
            cout << "    ||dx||/||x||/dt = " << setw(13) << twoNorm << endl;
        }
         
    }
    logger.cleanup();

    delete solver;
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
    else if ( modelName == "linearperiodic" ) {
        type = LINEARPERIODIC;
    }
    else if ( modelName == "sfd" ) {
        type = SFD;
    }
    else {
        cerr << "Unrecognized model: " << modelName << endl;
        type = INVALID;
    }
    return type;
}

Scheme::SchemeType str2scheme( string schemeName ) {
    Scheme::SchemeType type;
    MakeLowercase( schemeName );
    if ( schemeName == "euler" ) {
        type = Scheme::EULER;
    }
    else if ( schemeName == "ab2" ) {
        type = Scheme::AB2;
    }
    else if ( schemeName == "rk3" ) {
        type = Scheme::RK3;
    }
    else if ( schemeName == "rk3b" ) {
        type = Scheme::RK3b;
    }
    else {
        cerr << "Unrecognized integration scheme: " << schemeName;
        cerr << "    Exiting program." << endl;
        exit(1);
    }
    return type;
}


