// header file for storing integration schemes
// try to fix the miisalignment
#include "IBSolver.h"

#include "Geometry.h"
#include "ProjectionSolver.h"
#include "ConjugateGradientSolver.h"
#include "CholeskySolver.h"
#include "BoundaryVector.h"
#include "Grid.h"
#include "State.h"
#include "TangentSE2.h"
#include "VectorOperations.h"
#include "StateVector.h"

#include <iostream>
#include <fstream>

#include <math.h>
#include <string>


using std::ostream;
using std::ofstream;

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
		 cout << " updatesolver_NSsolver ... " << endl;
	}
	cout << " IBPMrk time " + std::to_string(i) + " ... " << x.time + _scheme.cn(i) * _dt << endl;
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

   Scalar AdjointIBSolver2::N(const State& x) const {
      // Get appropriate phase of periodic base flow
      //int k = _period - x.timestep - 1; // if not periodic
      int k = (_period - (x.timestep % _period)) % _period; // if periodic
			double t1 = remainder(x.time,((double)_period * _dt));
			if (t1 < 0){
				t1 = (double)_period*_dt+t1;
			}
		  double time = ((double)_period * _dt - t1); // if periodic
	  // double t2 = ((double)_period * _dt - t1); // if periodic
		// double time = remainder(t2,((double)_period * _dt));
		// if (time < 0){
		// 	time = (double)_period*_dt+time;
		// }
       cout << "At time step " << x.timestep << ", xtime = " << x.time << ", phase k = " << k << ", time = " << time << endl;
      // cout << "t1 = " << t1 << endl;
//////////
double xbod, ybod, theta; // position of rigid body
double xdot, ydot, thetadot; // velocity of rigid body
double xC = 0.0;
double yC = 0.0;
TangentSE2 gg = _motion->getTransformation(time);
gg.getPosition(xbod, ybod, theta);
gg.getVelocity(xdot, ydot, thetadot);
//
// qis stored for each (integer) time step from the forward simulation, but for the adjoint backward propergation using RK3,
// we need to have q at intermediate steps, so we use interpolation here.
Flux q(_x0periodic[k].q);
Flux q_plus(_x0periodic[k].q);
if (k < 1){
	Flux q_plus(_x0periodic[_period-1].q);
}
else{
	Flux q_plus(_x0periodic[k-1].q);   /// k+1
}

Flux q_temp(_x0periodic[k].q);

q_temp = q + (q_plus - q)/_dt * (x.time - x.timestep * _dt);

Flux q_pot(_x0periodic[k].q);
double dx = _grid.Dx();
Flux::index ind;
for (int lev=0; lev<q_pot.Ngrid(); ++lev) {

		// double xdot = -sin(theta);
		// double ydot = -cos(theta);
		// double xdot = 0;
		// double ydot = 1;
		double xdiff, ydiff;


		// double dx = dqda.Dx(lev);
		for(ind = q_pot.begin(X); ind != q_pot.end(X); ++ind) {
				xdiff = q_pot.x(lev,ind) - xC;
				ydiff = q_pot.y(lev,ind) - yC;
				// dqda(lev,ind) = (xdot -thetadot*ydiff)*dx;
				// q_pot(lev,ind) = (cos(theta) + thetadot * ydiff - ydot * sin(theta)) * dx;
				q_pot(lev,ind) = (thetadot * ydiff) * dx * 0;
		}
		for(ind = q_pot.begin(Y); ind != q_pot.end(Y); ++ind) {
				xdiff = q_pot.x(lev,ind) - xC;
				ydiff = q_pot.y(lev,ind) - yC;
				// dqda(lev,ind) = (ydot + thetadot*xdiff)*dx;
				// q_pot(lev,ind) = (-sin(theta) - thetadot * xdiff - ydot * cos(theta)) * dx;
				q_pot(lev,ind) = (-thetadot * xdiff) * dx * 0;
		}
}

// omega is stored for each (integer) time step from the forward simulation, but for the adjoint backward propergation using RK3,
// we need to have omega at intermediate steps, so we use interpolation here.
Flux q0_temp = q_temp + q_pot;
State x0_temp = _x0periodic[k];
State x0 = _x0periodic[k];
State x0_plus = _x0periodic[k];
if (k < 1){
	x0_plus = _x0periodic[_period-1];
}
else{
	x0_plus = _x0periodic[k-1];
}
x0_temp.omega = x0.omega + (x0_plus.omega - x0.omega)/_dt * (x.time - x.timestep * _dt);
// Scalar omega_tem = -(-2 * thetadot - _x0periodic[k].omega);
// Scalar omega_tem(_x0periodic[k].omega);
// omega_tem += 2;
// _x0periodic[k].omega = -(-2 * thetadot - _x0periodic[k].omega);
// _x0periodic[k].omega += _x0periodic[k].omega;
// cout << " thetadot1_solver ... " << thetadot << endl;
/////////

			// Scalar g = Laplacian( CrossProduct( _x0periodic[k].q, x.q )); // 1/delta^2   L = -C^TC   + q_pot !!!!!!!
			Scalar g = Laplacian( CrossProduct( q0_temp, x.q )); // 1/delta^2   L = -C^TC   + q_pot !!!!!!!
      g -= Curl( CrossProduct( x.q, x0_temp.omega));  // C^T q = csi  C(C^TC)^(-1) w = q ////////////////////////
			// g -= Curl( CrossProduct( x.q, omega_tem));  // C^T q = csi  C(C^TC)^(-1) w = q
      return g;
      // There is a term with dI/domega that should be added.
   }

   void AdjointIBSolver2::advanceSubstep( State& x, const Scalar& nonlinear, int i ) {
      // If the body is moving, update its position
      if ( _model.isTimeDependent() ) {
	 _model.updateOperators( x.time + _scheme.cn(i) * _dt ); /// need to be changed to _model.updateOperators( time - _scheme.cn(i) * _dt ); for moving body
	  cout << " updatesolver_ADJsolver ... " << endl;
      }
			 cout << " time_solver ... " << x.time + _scheme.cn(i) * _dt << endl;

			// _model.updateOperators( x.time + _scheme.cn(i) * _dt ); //////////////////

      // Evaluate Right-Hand-Side (a) for first equation of ProjectionSolver
      Scalar a = Laplacian( x.omega );
      a *= 0.5 * _model.getAlpha() * ( _scheme.an(i) + _scheme.bn(i) );
      a += _scheme.an(i)*nonlinear;

			double alpha = _model.getAlpha();
			// cout << " i ... " << i <<endl;
			 cout << " alpha_solver ... " << alpha <<endl;
       cout << " _scheme.an ... " << _scheme.an(i) << endl;
			// cout << " _scheme.bn ... " << _scheme.bn(i) << endl;
			// cout << " _scheme.cn ... " << _scheme.cn(i) << endl;
			// cout << " dt ... " << _dt << endl;
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
      double xbod, ybod, theta; // position of rigid body
      double xdot, ydot, thetadot; // velocity of rigid body
      //int k = _period - x.timestep - 1; // if not periodic
      int k = (_period - (x.timestep % _period)) % _period; // if periodic
	  // double time = ((double)_period * _dt - (x.time % ((double)_period * _dt))) % ((double)_period * _dt); // if periodic
		double t1 = remainder(x.time,((double)_period * _dt));
		if (t1 < 0){
			t1 = (double)_period*_dt+t1;
		}
	  double time = ((double)_period * _dt - t1); // if periodic
		// double time = remainder(t2,((double)_period * _dt));
		// if (time < 0){
		// 	time = (double)_period*_dt+time;
		// }

		/////////////////////////////
		TangentSE2 g = _motion->getTransformation(time);
				// g = _motion->getTransformation(time + (1.0 + 1.0/3.0) * _dt);
				// cout << " rk time (0) ... " << time + (1.0 + 1.0/3.0) * _dt << endl;
				if (time < _dt/2){
					g = _motion->getTransformation((double)_period * _dt - _scheme.cn(i) * _dt);
					// cout << " rk time " + std::to_string(i) + " ... " << (double)_period * _dt - _scheme.cn(i) * _dt << endl;
					cout << " rk time " + std::to_string(i) + " ... " << (double)_period * _dt - _scheme.cn(i) * _dt << endl;
				}
				else{
					if (time - _scheme.cn(i) * _dt < 0.0){
						g = _motion->getTransformation(0.0);
						cout << " rk time " + std::to_string(i) + " ... " << 0.0 << endl;
					}
					else{
						g = _motion->getTransformation(time - _scheme.cn(i) * _dt);
						cout << " rk time " + std::to_string(i) + " ... " << time - _scheme.cn(i) * _dt << endl;
					}
					 
				}
			/////////////////////////

			// if (time < _dt/2){
			// 	g = _motion->getTransformation((double)_period * _dt - _scheme.cn(i) * _dt);
			// 	cout << " rk time " + std::to_string(i) + " ... " << (double)_period * _dt - _scheme.cn(i) * _dt << endl;
			// }
			// else{
			// 	g = _motion->getTransformation(time - _scheme.cn(i) * _dt);
			// 	cout << " rk time " + std::to_string(i) + " ... " << time - _scheme.cn(i) * _dt << endl;
			// }
		/////////////////////////////
			// TangentSE2 g = _motion->getTransformation(time);
      // if(i == 0){
			// 	if (time+_dt > _period * _dt){
			// 		g = _motion->getTransformation(_period * _dt - 1.0/3.0 * _dt);  // check the sign
			// 		cout << " rk time (0) ... " << _period * _dt - 1.0/3.0 * _dt << endl;
			// 	}
			// 	else {
			// 		// g = _motion->getTransformation(time + (1.0 + 1.0/3.0) * _dt);
			// 		// cout << " rk time (0) ... " << time + (1.0 + 1.0/3.0) * _dt << endl;
			// 		g = _motion->getTransformation(time + (1.0 - 1.0/3.0) * _dt);
			// 		cout << " rk time (0) ... " << time + (1.0 - 1.0/3.0) * _dt << endl;
			// 	}
			//
			// }
			// else if (i == 1){
			// 	if (time+_dt > _period * _dt){
			// 		g = _motion->getTransformation(_period * _dt - 3.0/4.0 * _dt);
			// 		cout << " rk time (1) ... " << _period * _dt - 3.0/4.0 * _dt << endl;
			// 	}
			// 	else{
			// 		g = _motion->getTransformation(time + (1.0 - 3.0/4.0) * _dt);
			//   	 cout << " rk time (1) ... " << time + (1.0 - 3.0/4.0) * _dt << endl;
			// 	}
			// }
			// else if (i == 2){
			// 	g = _motion->getTransformation(time);
			// 	cout << " rk time (2) ... " << time << endl;
			// }
      /////////////////////////
			/////////////////////////////
				// TangentSE2 g = _motion->getTransformation(time);
				// if(i == 0){
				// 	if (time+_dt > _period * _dt){
				// 		g = _motion->getTransformation(0 * _dt + (1.0 - 1.0/3.0) * _dt);  // check the sign
				// 		cout << " rk time (0) ... " << 0 * _dt + (1.0 - 1.0/3.0) * _dt << endl;
				// 	}
				// 	else {
				// 		// g = _motion->getTransformation(time + (1.0 + 1.0/3.0) * _dt);
				// 		// cout << " rk time (0) ... " << time + (1.0 + 1.0/3.0) * _dt << endl;
				// 		g = _motion->getTransformation(time + (1.0 - 1.0/3.0) * _dt);
				// 		cout << " rk time (0) ... " << time + (1.0 - 1.0/3.0) * _dt << endl;
				// 	}
				//
				// }
				// else if (i == 1){
				// 	if (time+_dt > _period * _dt){
				// 		g = _motion->getTransformation(0 * _dt + (1.0 - 3.0/4.0) * _dt);
				// 		cout << " rk time (1) ... " << 0 * _dt + (1.0 - 3.0/4.0) * _dt << endl;
				// 	}
				// 	else{
				// 		g = _motion->getTransformation(time + (1.0 - 3.0/4.0) * _dt);
				// 		 cout << " rk time (1) ... " << time + (1.0 - 3.0/4.0) * _dt << endl;
				// 	}
				// }
				// else if (i == 2){
				// 	g = _motion->getTransformation(time);
				// 	cout << " rk time (2) ... " << time << endl;
				// }
				/////////////////////////
      g.getPosition(xbod, ybod, theta);
      g.getVelocity(xdot, ydot, thetadot);

			// cout << " xbod_solver ... " << xbod << endl;
			// cout << " ybod_solver ... " << ybod << endl;
			// cout << " xdot_solver ... " << xdot << endl;
			// cout << " ydot_solver ... " << ydot << endl;
			// cout << " thetadot_solver ... " << thetadot << endl;
			// cout << " theta_solver ... " << theta << endl;


			////////////////////// only for test
			// if (i == 0){
			// 	// read theta
			// 	std::ifstream file_theta0("../examples/theta_solver0.txt");
			// 	double theta_temp0[(_period + 1)];
			// 	// Update Error
			// 	for (int count = 0; count <= _period; ++count)
			// 		{
			// 			file_theta0 >> theta_temp0[count];
			// 		}
			// 	// Save theta
			// 	theta_temp0[k] = theta;
			// 	ofstream myfile_theta0("../examples/theta_solver0.txt");
			// 	if (myfile_theta0.is_open())
			// 	{
			// 		// for (int count = 0; count < size; count++) {
			// 		for (int count = 0; count <= _period; ++count) {
			// 			myfile_theta0 << theta_temp0[count] << " ";
			// 		}
			// 		myfile_theta0.close();
			// 	}
			// 	else cout << "Unable to open file";
			// 	//
			// 	// read y
			// 	std::ifstream file_y0("../examples/y_solver0.txt");
			// 	double y_temp0[(_period + 1)];
			// 	// Update Error
			// 	for (int count = 0; count <= _period; ++count)
			// 		{
			// 			file_y0 >> y_temp0[count];
			// 		}
			// 	// Save y
			// 	y_temp0[k] = ybod;
			// 	ofstream myfile_y0("../examples/y_solver0.txt");
			// 	if (myfile_y0.is_open())
			// 	{
			// 		// for (int count = 0; count < size; count++) {
			// 		for (int count = 0; count <=_period; ++count) {
			// 			myfile_y0 << y_temp0[count] << " ";
			// 		}
			// 		myfile_y0.close();
			// 	}
			// 	else cout << "Unable to open file";
			//
			// }
			//
			// else if (i == 1){
			// 	// read theta
			// 	std::ifstream file_theta1("../examples/theta_solver1.txt");
			// 	double theta_temp1[(_period + 1)];
			// 	// Update Error
			// 	for (int count = 0; count <= _period; ++count)
			// 		{
			// 			file_theta1 >> theta_temp1[count];
			// 		}
			// 	// Save theta
			// 	theta_temp1[k] = theta;
			// 	ofstream myfile_theta1("../examples/theta_solver1.txt");
			// 	if (myfile_theta1.is_open())
			// 	{
			// 		// for (int count = 0; count < size; count++) {
			// 		for (int count = 0; count <= _period; ++count) {
			// 			myfile_theta1 << theta_temp1[count] << " ";
			// 		}
			// 		myfile_theta1.close();
			// 	}
			// 	else cout << "Unable to open file";
			// 	//
			// 	// read y
			// 	std::ifstream file_y1("../examples/y_solver1.txt");
			// 	double y_temp1[(_period + 1)];
			// 	// Update Error
			// 	for (int count = 0; count <= _period; ++count)
			// 		{
			// 			file_y1 >> y_temp1[count];
			// 		}
			// 	// Save y
			// 	y_temp1[k] = ybod;
			// 	ofstream myfile_y1("../examples/y_solver1.txt");
			// 	if (myfile_y1.is_open())
			// 	{
			// 		// for (int count = 0; count < size; count++) {
			// 		for (int count = 0; count <=_period; ++count) {
			// 			myfile_y1 << y_temp1[count] << " ";
			// 		}
			// 		myfile_y1.close();
			// 	}
			// 	else cout << "Unable to open file";
			//
			// }
			//
			//
			//
			// else if (i == 2){
			// 	// read theta
			// 	std::ifstream file_theta2("../examples/theta_solver2.txt");
			// 	double theta_temp2[(_period + 1)];
			// 	// Update Error
			// 	for (int count = 0; count <= _period; ++count)
			// 		{
			// 			file_theta2 >> theta_temp2[count];
			// 		}
			// 	// Save theta
			// 	theta_temp2[k] = theta;
			// 	ofstream myfile_theta2("../examples/theta_solver2.txt");
			// 	if (myfile_theta2.is_open())
			// 	{
			// 		// for (int count = 0; count < size; count++) {
			// 		for (int count = 0; count <= _period; ++count) {
			// 			myfile_theta2 << theta_temp2[count] << " ";
			// 		}
			// 		myfile_theta2.close();
			// 	}
			// 	else cout << "Unable to open file";
			// 	//
			// 	// read y
			// 	std::ifstream file_y2("../examples/y_solver2.txt");
			// 	double y_temp2[(_period + 1)];
			// 	// Update Error
			// 	for (int count = 0; count <= _period; ++count)
			// 		{
			// 			file_y2 >> y_temp2[count];
			// 		}
			// 	// Save y
			// 	y_temp2[k] = ybod;
			// 	ofstream myfile_y2("../examples/y_solver2.txt");
			// 	if (myfile_y2.is_open())
			// 	{
			// 		// for (int count = 0; count < size; count++) {
			// 		for (int count = 0; count <=_period; ++count) {
			// 			myfile_y2 << y_temp2[count] << " ";
			// 		}
			// 		myfile_y2.close();
			// 	}
			// 	else cout << "Unable to open file";

			// }
			////////////////////////////
			///////////////////////////

      int n = _model.getNumPoints();
      double dx2 = _grid.Dx() * _grid.Dx();
      BoundaryVector b = _model.getConstraints();
      b = 0;
			// cost means cost function, 0 is drag, 1 is effeciency
		  if (_cost == 0){
				for (int j = 0; j < n; ++j) {
					////////////
					// -dg/dfx -dg/dfy
					b(X,j) += -2 * cos(theta) * dx2 / ((double)_period*_dt);     // pitchplunge for min drag + somthing else  ////// change sign!!!!!!
					b(Y,j) += 2 * sin(theta) * dx2 / ((double)_period*_dt);
					////////////
				}
			}


      else if (_cost == 1){
				State x0_temp = _x0periodic[k];
				State x0 = _x0periodic[k];
				State x0_plus = _x0periodic[k];
				if (k < 1){
					x0_plus = _x0periodic[_period-1];
				}
				else{
					x0_plus = _x0periodic[k-1];
				}
				x0_temp.f = x0.f + (x0_plus.f - x0.f)/_dt * (_scheme.cn(i) * _dt);

				// define the location vector of the boundary
						BoundaryVector boundary = _model.getPoints();
						double F_x = 0.;    // define the integral of fx * dx2
						double F_y = 0.;
						double M_z = 0.;

						for (int j = 0; j < n; ++j){
							F_y = 2 * x0_temp.f(1,j) * dx2 + F_y;
							F_x = 2 * x0_temp.f(0,j) * dx2 + F_x;
							M_z = 2 * x0_temp.f(1,j) * (boundary(X,j) - 0.0) * dx2 + M_z; // 0.0 should be replaced by xC in future
						}

						// cout << " F_y ... " << F_y << endl;
						// cout << " F_x ... " << F_x << endl;
						// cout << " M_z ... " << M_z << endl;


						for (int j = 0; j < n; ++j) {
							////////////////////
              b(X,j) = -(1/_CP * 2 * cos(theta) * dx2 / ((double)_period * _dt) - _CT / (pow(_CP,2)) * (-2) * sin(theta) * dx2 * ydot / ((double)_period * _dt));
							b(Y,j) = -(1/_CP * (-2) * sin(theta) * dx2 / ((double)_period * _dt) - _CT / (pow(_CP,2)) * (-2 * cos(theta) * dx2 * ydot - 2 * boundary(X,j) * dx2 * thetadot) / ((double)_period * _dt));

              // // -dg/dfx -dg/dfy
							// b(X,j) = -(-2 * sin(theta) * dx2 * ydot + 2 * cos(theta) * dx2) / ((double)_period*_dt);  // /(period*dt)
							//
							// b(Y,j) = -(-2 * cos(theta) * dx2 * ydot - 2 * (boundary(X,j)-0.0) * dx2 * thetadot - 2 * sin(theta) * dx2) / ((double)_period*_dt);
							// //////////// *2
					   //  b(X,j) = -2 * cos(theta) * dx2 * (-cos(theta) * F_y * ydot - sin(theta) * F_x * ydot - M_z * thetadot)/(pow(-cos(theta) * F_y * ydot - sin(theta) * F_x * ydot - M_z * thetadot,2)) / _period\
					 	 //          -(-F_x * cos(theta) + F_y * sin(theta)) * (-sin(theta) * ydot * dx2 * 2)/(pow(-cos(theta) * F_y * ydot - sin(theta) * F_x * ydot - M_z * thetadot,2)) / _period;
						 //
					 	 //  b(Y,j) = 2 * sin(theta) * dx2 * (-cos(theta) * F_y * ydot - sin(theta) * F_x * ydot - M_z * thetadot)/(pow(-cos(theta) * F_y * ydot - sin(theta) * F_x * ydot - M_z * thetadot,2)) / _period\
					 	 //          -(-F_x * cos(theta) + F_y * sin(theta))*(-cos(theta) * ydot * 2 - (boundary(X,j)-0.0) * thetadot * 2) * dx2/(pow(-cos(theta) * F_y * ydot - sin(theta) * F_x * ydot - M_z * thetadot,2)) / _period;
					 	 // ////////////

						}
						// cout << " b(X,j) check ... " << (pow(-cos(theta) * F_y * ydot - sin(theta) * F_x * ydot - M_z * thetadot,2)) << endl;
						// cout << " cos(theta) ... " << cos(theta) << endl;
						// cout << " sin(theta) ... " << sin(theta) << endl;
						// cout << " ydot ... " << ydot << endl;
						// cout << " thetadot ... " << thetadot << endl;
						// cout << " _CT ... " << _CT << endl;
						// cout << " _CP ... " << _CP << endl;



			}

		  // cout << " dx2 ... " << dx2 << endl;
      // b = 0; // b = (dI/df)^T, need to program this somehow

      // Call the ProjectionSolver to determine the vorticity and forces
	 _solver[i]->solve( a, b, x.omega, x.f );

      // Update the state, for instance to compute the corresponding flux
	 _model.refreshState( x ); // this may add  some base flow to flux, which is not needed for adjoint. Check this.
      _Nprev = nonlinear;

      if( _oldSaved == false ) {
	 _oldSaved = true;
      }
   }

Scalar LinearizedPeriodicIBSolver::N(const State& x) const {
	int k = x.timestep % _period;
	// cout << "At time step " << x.timestep << ", phase k = " << k << endl;
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
