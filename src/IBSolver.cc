// header file for storing integration schemes
#include "IBSolver.h"

#include "Geometry.h"
#include "ProjectionSolver.h"
#include "ConjugateGradientSolver.h"
#include "CholeskySolver.h"
#include "BoundaryVector.h"
#include "Grid.h"
#include "State.h"
#include "VectorOperations.h"
#include <string>

namespace ibpm {

IBSolver::IBSolver( 
	const Grid& grid, 
	NavierStokesModel& model,
	double dt, 
	Scheme::SchemeType scheme
	) :
	_grid( grid ),
	_scheme( scheme ),
	_name( _scheme.name( ) ),
	_dt( dt ),
	_model( model ),
	_Nprev( grid ),
	_Ntemp( grid ), 
	_oldSaved( false ),
	_solver( _scheme.nsteps() ),
    _tol( 1e-7) {	
		createAllSolvers();
	}
	
IBSolver::IBSolver( 
    const Grid& grid, 
    NavierStokesModel& model,
    double dt, 
    Scheme::SchemeType scheme,
    double tol
    ) :
    _grid( grid ),
    _scheme( scheme ),
    _name( _scheme.name( ) ),
    _dt( dt ),
    _model( model ),
    _Nprev( grid ),
    _Ntemp( grid ), 
    _oldSaved( false ),
    _solver( _scheme.nsteps() ),
    _tol( tol ) {	
        createAllSolvers();
}
	
IBSolver::~IBSolver() {
	deleteAllSolvers();
}
	
string IBSolver::getName() {
	return _name;
}
    
double IBSolver::getTimestep() {
    return _dt;
}
	
void IBSolver::init() {	
	for ( int i = 0; i < _scheme.nsteps(); i++ ) {
		_solver[i] -> init();
	}
}
    
void IBSolver::reset() {	
    _oldSaved = false;
}    
	
bool IBSolver::load(const string& basename) {
	bool successInit = false;
	bool successTemp = true;
	for ( int i = 0; i < _scheme.nsteps(); i++ ) {
		char num[256];
		sprintf( num, "%02d", i+1 );
		string filename = basename + "_" + num;
		successTemp = _solver[i] -> load( filename ) && successTemp;
        if ( i == 0 ) { 
            successInit = true; 
        }
	}

	return successInit && successTemp;
}
	
bool IBSolver::save(const string& basename) {
	bool successInit = false;
	bool successTemp = true;
	for ( int i = 0; i < _scheme.nsteps(); i++ ) {
		char num[256];
		sprintf( num, "%02d", i+1 );
		string filename = basename + "_" + num;
		
		successTemp = _solver[i] -> save( filename ) && successTemp;
		if ( i == 0 ) { 
            successInit = true; 
        }
	}
	
	return successInit && successTemp;
}
	
void IBSolver::createAllSolvers() {
	for ( int i = 0; i < _scheme.nsteps(); i++ ) {
		_solver[i] = createSolver( ( _scheme.an(i) + _scheme.bn(i) )*_dt );
	}
}
	
void IBSolver::deleteAllSolvers() {
	for (unsigned int i = 0; i < _solver.size(); i++) {
		delete _solver[i];
	}
}
	
ProjectionSolver* IBSolver::createSolver(double beta) {
	// Check whether all bodies are stationary
	//      If so, return a CholeskySolver
	//      If not, return a ConjugateGradientSolver
	if ( _model.geTimeDependent() ) {
		cerr << "Using ConjugateGradient solver for projection step" << endl
		<< "  tolerance = " << _tol << endl;
		return new ConjugateGradientSolver( _grid, _model, beta, _tol );    
	}
	else {
		cerr << "Using Cholesky solver for projection step" << endl;
		return new CholeskySolver( _grid, _model, beta );
	}
}
    
void IBSolver::setTol( double tol ) {
    _tol = tol;
    createAllSolvers();
}

void IBSolver::advance( State& x ) {	
	for ( int i = 0; i < _scheme.nsteps(); i++ ) {
		Scalar nonlinear = N(x); 
		advanceSubstep( x, nonlinear, i );
	}
    
	x.time += _dt;
	++x.timestep;
}      
	
void IBSolver::advance( State& x, const Scalar& Bu ) {
	for ( int i = 0; i < _scheme.nsteps(); i++ ) {
		Scalar nonlinear = N(x) + Bu; 
		advanceSubstep( x, nonlinear, i );
	}
    
    x.time += _dt;
	++x.timestep;
}     
	
void IBSolver::advanceSubstep( State& x, const Scalar& nonlinear, int i ) {    
	// If the body is moving, update the positions of the bodies
	if ( _model.isTimeDependent() ) {
		_model.updateOperators( x.time + _scheme.cn(i) * _dt );
	}		
    
	// Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
	Scalar a = Laplacian( x.omega );
	a *= 0.5 * _model.getAlpha() * ( _scheme.an(i) + _scheme.bn(i) );
	a += _scheme.an(i)*nonlinear;
	
	if ( _scheme.bn(i) != 0 ) {        
        // for ab2
		if ( _oldSaved == false ) {
			_Nprev = nonlinear;
		}
        
		a += _scheme.bn(i) * _Nprev;
	}
	
	a *= _dt;
	a += x.omega;

	// Evaluate Right-Hand-Side (b) for second equation of ProjectionSolver
	BoundaryVector b = _model.getConstraints();
    
	// Call the ProjectionSolver to determine the vorticity and forces
	_solver[i]->solve( a, b, x.omega, x.f );

	// Update the state, for instance to compute the corresponding flux
	_model.refreshState( x );	
	_Nprev = nonlinear;
    
    if( _oldSaved == false ) {
        _oldSaved = true;       
    }
}	
	
	
// ===================== //
// Derived class methods //
// ===================== //
	
Scalar NonlinearIBSolver::N(const State& x) const {
	Flux v = CrossProduct( x.q, x.omega );
	Scalar g = Curl( v );
	return g;
}
	
Scalar LinearizedIBSolver::N(const State& x) const {
	Flux v = CrossProduct( _x0.q, x.omega );
	v += CrossProduct( x.q, _x0.omega );
	Scalar g = Curl( v );
	return g;
}	
	
Scalar AdjointIBSolver::N(const State& x) const {
    Scalar g = Laplacian( CrossProduct( _x0.q, x.q ));
	g -= Curl(CrossProduct( x.q, _x0.omega ));
	return g;
}	
	
Scalar LinearizedPeriodicIBSolver::N(const State& x) const {
	int k = x.timestep % _period;
	cout << "At time step " << x.timestep << ", phase k = " << k << endl; 
	Flux v = CrossProduct( _x0periodic[k].q, x.omega );
	v += CrossProduct( x.q, _x0periodic[k].omega );
	Scalar g = Curl( v );	
	return g;		
}
	

// =========== //
// SFD methods //
// =========== //
	
Scalar SFDSolver::N(const State& x) const {
	Flux v = CrossProduct( x.q, x.omega );
	Scalar g = Curl( v );
	Scalar temp( x.omega );  // because x is const here...hmmm
	g -= _chi * ( temp - _xhat.omega );
	return g;
}
	

void SFDSolver::advanceSubstep( State& x, const Scalar& nonlinear, int i ) {
    assert( x.time == _xhat.time );
    
	// Initialize _xhat if necessary, save current vorticity field
	if ( _xhatSaved == false ) {
		_xhat = x;
		_xhatSaved = true;
	}
	_omegaTemp = x.omega;
	
	// Advance state x
	IBSolver::advanceSubstep( x, nonlinear, i );
	
	// Advance state _xhat
	Scalar rhs = ( _omegaTemp - _xhat.omega ) / _Delta;
	Scalar a = _scheme.an(i)*rhs;
	
	if ( _scheme.bn(i) != 0 ) {
		if ( _rhsSaved == false ) {
			_rhsPrev = rhs;
		}
		
		a += _scheme.bn(i) * _rhsPrev;
	}
	
	a *= _dt;
	_xhat.omega += a;
    
    if ( i == _scheme.nsteps()-1 ) {
        _xhat.time += _dt;
        _xhat.timestep++;
    }
    
	_rhsPrev = rhs;
    
    if( _rhsSaved == false ) {
        _rhsSaved = true;       
    }
}
 
void SFDSolver::saveFilteredState( string outdir, string name, string numDigitInFileName ) { 
    string formatString = outdir+name+numDigitInFileName+".bin"+"_xhat";
    char filename[256];
    sprintf( filename, formatString.c_str(), _xhat.timestep );
    _xhat.save( filename );
}
    
void SFDSolver::loadFilteredState( string icFile ) { 
    string xhatFile = icFile+"_xhat";
    _xhat.omega = 0.;
    _xhat.f = 0.;
    _xhat.q = 0.;
	if (xhatFile != "_xhat") {	
	    cout << "Loading initial condition from file: " << xhatFile << endl;
        if ( ! _xhat.load(xhatFile) ) {
            cout << "  (failed: setting xhat = x)" << endl;
        }
    }
    else {
        cout << "Setting xhat = x" << endl;
    }
    
    _xhatSaved = true;
}

} // ibpm
