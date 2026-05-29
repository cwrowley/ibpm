// A working version of ADJ optimization, using newton-gmres on the adj part,
// rk3 is acceptable here
// the first point of theta and y can be updated correctly and the limit cycle is solved correclty
// This is a working code without constraint

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
#include <Eigen>
#include <Dense>
#include <Core>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>
#include <exception>
#include <set>
#include <chrono>
#include <ctime>
#include <ratio>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
//#include <iterator>
//#include <limits>
//#include <vector>
//#include <iterator>
//#include <LBFGS.h>
#include <LBFGSB.h>

#include "StateVector.h"
#include "Gmres.h"
#include "Newton.h"
#include "NewtonArmijo.h"
#include "exceptions.h"
#include "mpi.h"
#include "MotionGenerator.h"


//using namespace std;
using std::ostream;
using std::cout;
using std::string;
using std::endl;
using std::cerr;
using std::setw;
using namespace std::chrono;
using namespace ibpm;
using std::ofstream;
using std::sort;
//using namespace Eigen;
//using namespace newton_ns; //

using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
using Scal = ibpm::Scalar; ////////////////////////// Scalar

using namespace LBFGSpp;


// Define global variable parameter
int size1 = 800;     ////period/dt   400
int size2 = 2;    // 2 pitch and heave
int size3 = 1;    // 2  number of Fourier modes
double Fourier_mode = 0;    // 0 for y = sum (An * sin(t + Bn)), theta = sum (Cn * sin(t + Dn)); 1 for y = sum (An * sin + Bn * cos), theta = sum (Cn * sin + Dn * cos)
int cost = 1; // costfunction, cost = 0 is for drag, cost = 1 is for efficiency
int num_cycles = 1; // number of cycles
int Ini_step = 8; //initial cycles for newton solver
int k_global;
int l_global;
const int num_iter = 20;
const int num_iter_linsearch = 20;
double drag_global[num_iter];
double alpha_step = 0.09; // learning rate, which is not used here

string outdir = "Motion_Parameters_eff_BFGS_dt_1mode_Re100_MEM";
string outdir_PBF;
string outdir_PBF_global = "PBF_eff_BFGS_dt_1mode_Re100_MEM";
string outdir_ADJ = "ADJ_eff_BFGS_dt_1mode_Re100_MEM";

///////
// optimization parameters
double theta_up = 1.0;
double A_lim = 2.0;  // 1.0 heave
double B_lim = 3.14;  // 3.14 heave
double C_lim = 1.0;  // 0.4 pitch
double D_lim = 3.14;  // 3.14 pitch
double dt_lim_up = 0.02;
double dt_lim_low = 0.0045;


Vec A_0(size3);
Vec B_0(size3);
Vec C_0(size3);
Vec D_0(size3);
Vec z_m_save(size1+1);
Vec L_save(size1+1);
Vec D_save(size1+1);


void initialization (){
	for (int i = 0; i<size3; ++i){
		A_0(i) = 1.2;    // 0.4 add a phase differnce otherwise denominater could be 0 0.7
		B_0(i) = 1.45837;    // 1.57, 0.1 1.57
		C_0(i) = 0.111184;     // 0.1  0.4
		D_0(i) = 0.00690651;      // 0.  0.0
		/////////////////
		// Turn this off when validating the ADJ code against FD
		 if (i > 0){
		 	A_0(i) = 0.;
		 	B_0(i) = 0.;
		 	C_0(i) = 0.;
		 	D_0(i) = 0.;
		 }
		/////////////////
	}
}


/////////

Vec A(size3);
Vec B(size3);
Vec C(size3);
Vec D(size3);

Vec dA(size3);
Vec dB(size3);
Vec dC(size3);
Vec dD(size3);
double dT;  // gradient with respect to period
double dT_nonperiodic;
double g_T; // the value of cost function at the time instant T, T is period
double G_T;
double h_T;

Mat MotionParameter(size1 * Ini_step + 1,size2);
Vec MotionTime(size1 * Ini_step + 1);
Mat MotionParameterTime(size1 * Ini_step + 1,size2+1);

Mat MotionParameterTPerturb(size1 * Ini_step + 1,size2);
Vec MotionTimeTPerturb(size1 * Ini_step + 1);
Mat MotionParameterTimeTPerturb(size1 * Ini_step + 1,size2+1);

//
double meandrag[num_iter * num_iter_linsearch];
double drag_int = 0;
double dt = 0.01;    // 0.01
// define dt
double delta_T = dt / (double) size1;

// double Newton_tol = 1.e-3;
double Reynolds = 100;
double magnitude = 1;
int nx = 100;   // 100
int ny = 100;   // 100
int ngrid = 4;
bool TecplotAllGrids = 0;
double length = 2.0; // 4.0
double xOffset = -0.5; // -1
double yOffset = -1.0; // -2
double dy = abs(yOffset) * 2 / ny;
double xC = 0.0;
double yC = 0.0;
double CT = 0.0;
double CP = 0.0;
int save_flag = 0; // don't save the flow for periodic orbit solver
int T_perturb = 0; // create new folder for T_perturb so if won't conflict with the ADJ solver when parallelized

///////////////////////////////////////
// define vector<State> to store the base flow
vector<State> x0_glob;
////////////////////////////////

const double PI = 4. * atan(1.);
double twopi = 8. * atan(1.);


// Define innerproduct
inline double innerproduct(const Vec& x, const Vec& y) {
	return y.dot(x);
}
// Define norm
inline double norm(const Vec& x) {
	return x.norm();
}
// Define mean
inline double getAverage(const Vec& arr, int size) {
  int i = 0;
	double sum = 0;
  double avg;

   for (i = 0; i < size; ++i) {
      sum += arr(i);
   }
   avg = double(sum) / size;

   return avg;
}


//
// Function to compute pitching for one cycle
inline Vec F(const Vec& Omega_in) {
	cout << "Pitching flat plate example\n";

	// lift and drag integral of drag
	double lift = 0.;
	double lift_pre = 0;
	double drag = 0.;
	drag_int = 0;
	CT = 0;
	CP = 0;
	z_m_save.resize(size1 + 1);
	L_save.resize(size1 + 1);
	D_save.resize(size1 + 1);
	// Setup grid
	Grid grid(nx, ny, ngrid, length, xOffset, yOffset);
  //time step

	// Load period for phase condition
	std::ifstream input_fileP("../examples/period.txt");
	double period;
	input_fileP >> period;

  period = size1; // remove if period not known

	// Load PhaseCal_flag
	std::ifstream input_PF("../examples/PhaseCal_flag.txt");
	int PhaseCal_flag;
	input_PF >> PhaseCal_flag;


	// Make a flat plate, length 1, with center at 1/4 chord
	RigidBody plate;
	plate.addLine(0, 0, 1, 0, grid.Dx());
	plate.setCenter(0, 0);

	// Set the motion to pitching, amplitude = 0.25, period 10 time units


	double pitch_in[size1 + 1];
	double plunge_in[size1 + 1];
	double motiontime_in[size1 + 1];

if(T_perturb == 0){
	for(int i = 0; i<size1 + 1; ++i){
	// for(int i = 0; i<size1+1; ++i){
				pitch_in[i] =  MotionParameter(i,0);
				plunge_in[i] =  MotionParameter(i,1);
				motiontime_in[i] = MotionTime(i);
	}
}
else{
	for(int i = 0; i<size1 + 1; ++i){
	// for(int i = 0; i<size1+1; ++i){
				pitch_in[i] =  MotionParameterTPerturb(i,0);
				plunge_in[i] =  MotionParameterTPerturb(i,1);
				motiontime_in[i] = MotionTimeTPerturb(i);
	}
}

	// create dqdb
  std::ifstream file_tc("../examples/theta_check.txt");
  double t_t[size1+1];
  // Update Error
  for (int i = 0; i < size1+1; ++i)
 	 {
 		 file_tc >> t_t[i];
 	 }
  // Save dqdb
  ofstream myfile_thetacheck("../examples/theta_check.txt");
  if (myfile_thetacheck.is_open())
  {
 	 // for (int count = 0; count < size; count++) {
 	 for (int count = 0; count < size1+1; ++count) {
 		 myfile_thetacheck << pitch_in[count] << " ";
 	 }
 	 myfile_thetacheck.close();
  }
  else cout << "Unable to open file";
	////////////////////

	//double freq = 1;
    //bool ubf = true;
	int TT = size1;
	int end_index = TT;
	PitchPlungeT_mod motion(
		dt,
		motiontime_in,
		pitch_in,
		plunge_in,
		end_index
	);
	plate.setMotion(motion);
	Geometry geom;
	geom.addBody(plate);
	//geom.moveBodies(0);

	// Setup equations to solve

	double alpha = 0;  // angle of background flow
	BaseFlow q_potential(grid, magnitude, alpha);

/////////////
	Motion* m = geom.transferMotion(); // pull motion from first RigidBody object
	geom.transferCenter(xC, yC);  // pull center of motion from RigidBody object
	q_potential.setMotion(*m);
	q_potential.setCenter(xC, yC);
////////////
	cout << "Setting up Navier Stokes model..." << flush;
	NavierStokesModel model(grid, geom, Reynolds, q_potential);
	model.init();
	cout << "done" << endl;
	//cout << "amplitude = " << A << endl;
	//cout << "amplitude = " << amplitude << endl;

	// Setup timestepper
	// NonlinearIBSolver solver(grid, model, dt, Scheme::EULER);
	// NonlinearIBSolver solver(grid, model, dt, Scheme::AB2);
	NonlinearIBSolver solver(grid, model, dt, Scheme::RK3);
	solver.init();

	// Build the state variable, zero initial conditions
	State x(grid, geom.getNumPoints());
	x.omega = 0.;
	x.q = 0.;
	x.f = 0.;


	// Build a state structure to store the points coordinate on the body
	State Boundary(grid, geom.getNumPoints());
	Boundary.f = 0.1;
	// Convert Vec variables to State
	int N = 0;

	for (int k = 0; k<ngrid;++k){
		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				x.omega(k, i, j) = Omega_in(N);   //x.omega(0,i,j)
				N = N + 1;
			}

		}

		for (int i = 0; i <= nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				x.q(k, 0, i, j) = Omega_in(N);   //x.omega(0,i,j)
				N = N + 1;
			}

		}

		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j <= ny; ++j) {
				x.q(k, 1, i, j) = Omega_in(N);   //x.omega(0,i,j)
				N = N + 1;
			}

		}

	}


	for (int i = 0; i < geom.getNumPoints(); ++i) {
		x.f(0, i) = Omega_in(N);   //x.omega(0,i,j)
		N = N + 1;
	}

	for (int i = 0; i < geom.getNumPoints(); ++i) {
		x.f(1, i) = Omega_in(N);   //x.omega(0,i,j)
		N = N + 1;
	}


Logger logger;


// Create output directory, if does not already exist
mkdir( outdir_PBF.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

// Setup output routines
// if (save_flag == 1){
	OutputForce force(outdir_PBF + "PBF_force.dat");
	OutputTecplot tecplot(outdir_PBF + "PBF%03d.plt", "Pitching plate, step %03d", 1);
	if(TecplotAllGrids) tecplot.setFilename( outdir_PBF + "PBF%03d_g%01d.plt" );
	 // }

	 OutputRestart restart(outdir_PBF + "pitch%03d.bin");

	 //Output Tecplot file every few timesteps
	 // logger.addOutput(&tecplot, 1);
	 // logger.addOutput(&force, 1);

	 // save flow for one cycle
	 if (save_flag == 1){
		 x.load(outdir_PBF + "pitch000.bin");
		 logger.addOutput(&tecplot, 10 );
		 logger.addOutput(&force, 10 );
		 logger.addOutput(&restart, 1);
	 }
	 else{
	 // when runing for the periodic orbit, don't save the flow to accelerate te computation
		 logger.addOutput(&restart, 1000);
	 }

	 //logger.doOutput(x);

  logger.init();

	double zm_ini = 0;
	double theta_ini = motion.gettheta(0 * dt);

	// logger.doOutput(theta_ini, q_potential, x, zm_ini);
	logger.doOutput(x);


	// Step
	const double PI = 4. * atan(1.);
	//
	// initialize phase condition parameters
	double indicator;
	int StopCondition = 0;
	int I_ini = 0;
	int I_final = 0;
/*	double cycle = 1 / freq / dt; */
	int T = size1;
	int numSteps = (int)T;
	if (PhaseCal_flag == 0) {
	    numSteps = (int)T;    //+1?
	}
	else {
		numSteps = (int)T;
	}

	z_m_save.resize(size1+1);
	L_save.resize(size1+1);
	D_save.resize(size1+1);
	//cout << "numSteps: " << period << endl;
	//cout << "PhaseCal_flag: " << PhaseCal_flag << endl;
	for (int i = 0; i < numSteps; ++i) {             /// +1 or -1
		//double theta = 1 * sin(2 * PI / T * x.time);
	/*	cout << "step " << setw(4) << i
			<< "  time = " << setw(5) << x.time;*/
			////

		double xF, yF;
		x.computeNetForce(xF, yF);
		////////////////////
		double ydot = motion.getydot(x.time);
		double y = motion.gety(x.time);
		double theta = motion.gettheta(x.time);
		double thetadot = motion.getthetadot(x.time);

		////////////////////
		//q_potential.setAlphaMag(x.time);
		// alpha = q_potential.getAlpha();
		drag = xF * cos(theta) - yF * sin(theta);
		lift = xF * sin(theta) + yF * cos(theta);



					/////////////////////////////
		Boundary.f = plate.getPoints();

		double zM = 0.0;
		double drag_check = 0.0;
		for (int j = 0;j < geom.getNumPoints();++j) {
				zM = x.f(1, j) * grid.Dx() * grid.Dx() * (Boundary.f(0, j) - xC) + zM;
				drag_check = x.f(0, j) * grid.Dx() * grid.Dx() + drag_check;
		}
		z_m_save(i) = zM;
		L_save(i) = lift;
		D_save(i) = drag;
		cout << " step_aft " << setw(4) << i << " time_aft " << setw(5) << x.time << endl;
		cout << " lift_aft ..." << lift << " drag_aft ..." << drag << endl;
		// cout << " theta_aft ... " << theta << endl;
		// cout << " thetadot_aft ... " << thetadot << endl;
		// cout << " y_aft ... " << y << endl;
		// cout << " ydot_aft ... " << ydot << endl;
		////////////////////////////////
		// Compute drag integral
		CT = CT + drag * magnitude * 2 / ((double)period*dt) * dt;     //  /tau * 2
		CP = CP + (- lift * ydot - zM * thetadot) * 2 / ((double)period*dt) * dt;

		// CT = CT + drag * magnitude * 2 * dt;     //  /tau * 2
		// CP = CP + (- lift * ydot - zM * thetadot) * 2 * dt;

		g_T = 0.0;
		G_T = 0.0;
		h_T = 0.0;
    if (cost == 0){
			g_T = drag * magnitude / ((double)period*dt) * 2.;
		}

		else if (cost == 1){
			G_T = drag * magnitude / ((double)period*dt) * 2.;
			h_T = (- lift * ydot - zM * thetadot) / ((double)period*dt) * 2.;
			g_T = G_T / CP - CT/(pow(CP, 2.)) * (h_T);
			// g_T = G_T / (CP / ((double)period * dt)) - (CT / ((double)period * dt))/(pow(CP / ((double)period * dt),2)) * (h_T);
			// g_T = (drag * magnitude - lift * ydot - zM * thetadot) / ((double)period*dt) * 2;
		}
		// cout << " CT ... " << CT;
		// cout << " CP ... " << CP;
		// cout << " drag_int ... " << (CT / (double)numSteps) / (CP / (double)numSteps);
    // cout << " period ... " << period << endl;
		// cout << " ((double)period * dt) ... " << ((double)period * dt) << endl;
		// cout << " drag ... " << drag << endl;
		// cout << " G_T ... " << G_T << endl;
		// cout << " h_T ... " << h_T << endl;
		// cout << " g_T ... " << g_T << endl;
    // logger.doOutput(theta, q_potential, x, zM);
	solver.advance(x);
	  logger.doOutput(x);
    //////////////////////////////



    //logger.doOutput(theta, q_potential , x, zM)
    ////////////////////////////
		// Compute drag integral

		//cout << "CT..." << CT << endl;
		//cout << "CP..." << CP << endl;


		// if PhaseCal_flag=0, find the period
		if (PhaseCal_flag == 0) {
			if (lift * lift_pre < 0 && lift < lift_pre) {
				StopCondition = StopCondition + 1;
			}
			lift_pre = lift;

			if (StopCondition == 0) {
				I_ini = i;
			}
		/*	cout << "indicator: " << indicator << endl;
			cout << "StopCondition: " << StopCondition << endl;*/

		}

		//

	}

	// test only
	////////////////////////////////
	////////////////////////////////
		// save L
	ofstream myfile_L_check("../examples/" + outdir + "L_check.txt");
	if (myfile_L_check.is_open())
	{
		// for (int count = 0; count < size; count++) {
		for (int count = 0; count < numSteps; ++count) {
			myfile_L_check << L_save[count+1] << " ";
		}
		myfile_L_check.close();
	}
	else cout << "Unable to open file";

	// save D
	ofstream myfile_D_check("../examples/" + outdir + "D_check.txt");
	if (myfile_D_check.is_open())
	{
		// for (int count = 0; count < size; count++) {
		for (int count = 0; count < numSteps; ++count) {
			myfile_D_check << D_save[count + 1] << " ";
		}
		myfile_D_check.close();
	}
	else cout << "Unable to open file";

	
	///////////////////////////////
	//////////////////////////////
	///////////////////////////////
	// save period if PhaseCal_flag=0
	if (PhaseCal_flag == 0) {
		period = I_final - I_ini;
		period = size1; // remove if period not known
		ofstream file_PF;
		file_PF.open("../examples/period.txt");
		file_PF << period;
		file_PF.close();
		//cout << "period: " << period << endl;
	}
	//

	// Convert State variables in to Vec
	int M = 0;
	//Vec Omega(((nx-1) * (ny-1) * 3) * ngrid + geom.getNumPoints() * 2);
	Vec Omega(((nx - 1)* (ny - 1) + (nx + 1) * ny + nx * (ny + 1))* ngrid + geom.getNumPoints() * 2);

	for (int k = 0; k < ngrid; ++k){
		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				Omega(M) = x.omega(k, i, j);
				M = M + 1;
			}

		}

		for (int i = 0; i <= nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				Omega(M) = x.q(k, 0, i, j);
				M = M + 1;
			}

		}

		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j <= ny; ++j) {
				Omega(M) = x.q(k, 1, i, j);
				M = M + 1;
			}

		}

	}


	for (int i = 0; i < geom.getNumPoints(); ++i) {
		Omega(M) = x.f(0, i);
		M = M + 1;
	}

	for (int i = 0; i < geom.getNumPoints(); ++i) {
		Omega(M) = x.f(1, i);
		M = M + 1;
	}

	//for (int i = 0; i < _numPoints; ++i) {
	//	x.f(0, i) = Omega_in(N);   //x.omega(0,i,j)
	//	N = N + 1;
	//}

	//drag_int = drag_int / T;

	// drag_int = (CT / (double)numSteps) / (CP / (double)numSteps);

	logger.cleanup();
	// Save drag integral
	//ofstream file_drag_int;
	//file_drag_int.open("drag_int.txt");
	//file_drag_int << drag_int;
	//file_drag_int.close();
	//cout << "drag_in: " << drag_int << endl;


	// Compute the norm of x_final-x_initial
	double f_norm = norm(Omega - Omega_in);
	cout << "f_norm ... " << f_norm << endl;

	// Load step
	//std::ifstream input_file_step("../examples/step.txt");
	//int step;
	//input_file_step >> step;
	//// Load Error
	//std::ifstream file_Err("../examples/ErrOut.txt");
	//double Err[800];
	//int size = 800;
	//// Update Error
	//for (int i = 0; i < size; ++i)
	//	{
	//		file_Err >> Err[i];
	//	}

	//// Save Error
	//Err[step] = { f_norm };
	//ofstream myfile("../examples/ErrOut.txt");
	//if (myfile.is_open())
	//{
	//	// for (int count = 0; count < size; count++) {
	//	for (int count = 0; count < size; ++count) {
	//		myfile << Err[count] << " ";
	//	}
	//	myfile.close();
	//}
	//else cout << "Unable to open file";
	////cout << "Error" << Err[step] << endl;

	//step = step + 1;
	//// Save step
	//ofstream file;
	//file.open("step.txt");
	//file << step;
	//file.close();
	////cout << "Step_in: " << step << endl;

	return Omega - Omega_in;
	//return Omega;

}
// Newton and drag integral

inline double drag_int_fun(const Mat& parameter) {

	using namespace std::chrono;
	// Start timing
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	// Start Newton solver after a few IBPM cycles
	//
	using std::cout;
	using std::cerr;
	using std::endl;
	using std::string;
	using std::ofstream;

	// lift and drag
	double lift = 0.;
	double drag = 0.;
	// time step

	// Load all variables from parameter

	/*double plungePeriod_temp = round(plungePeriod / dt);
	plungePeriod = plungePeriod_temp * dt;*/
	////////////////////////////////////////////////
	// pitch and plunge
	int TT = (size1) * Ini_step;
	int end_index = TT;
	// double pitchPeriod_temp = round(TT / dt);
	// // TT = pitchPeriod_temp * dt;
  // int pitchPeriod = (int)pitchPeriod_temp;

	double pitch_in[size1 * Ini_step + 1];
	double plunge_in[size1 * Ini_step + 1];
	double motiontime_in[size1 * Ini_step + 1];


	for(int i = 0; i<size1 * Ini_step + 1; ++i){
				pitch_in[i] =  parameter(i,0);
				plunge_in[i] =  parameter(i,1);
				motiontime_in[i] = parameter(i,2);
	}

	// Setup grid


	Grid grid(nx, ny, ngrid, length, xOffset, yOffset);
	// Make a flat plate, length 1, with center at 1/4 chord
	RigidBody plate;
	plate.addLine(0, 0, 1, 0, grid.Dx());
	plate.setCenter(0, 0);

	// Set the motion to pitching, amplitude = 0.25, period 10 time units
	PitchPlungeT_mod motion(
		   dt,
			 motiontime_in,
	     pitch_in,
	     plunge_in,
	     end_index
	);
	// Save motion configuration
	//int size_para = sizeof(parameter_coef) / sizeof(parameter_coef[0]);


	/*const int size_p = 11;
	double para[11] = { amplitude0, Tp0, amplitude1, Tp1, amplitude2, Tp2, amplitude3, Tp3, amplitude4, Tp4, T };*/
	//cout << "T..." << T << endl;
	//cout << "Td..." << T << endl;
	//ofstream myfile_para;
	//myfile_para.open("../examples/para.txt");

	//for (int count = 0; count < size_p; count++) {
	//for (int count = 0; count < size_p; ++count) {
	//	myfile_para << para[count] << " ";
	//}
	//myfile_para.close();

	// Save step
	int step = Ini_step;
	ofstream file_step;
	file_step.open("../examples/step.txt");
	file_step << step;
	file_step.close();
	cout << "Step_in: " << step << endl;
	//Save Error
	const int size = 800;
	double Err[800] = {};

	ofstream myfile_ErrIc;
	myfile_ErrIc.open("../examples/ErrOut.txt");

	//for (int count = 0; count < size; count++) {
	for (int count = 0; count < size; ++count) {
		myfile_ErrIc << Err[count] << " ";
	}
	myfile_ErrIc.close();
	// Save drag integral
	//drag_int = 0;
	//ofstream file_drag_int;
	//file_drag_int.open("../examples/drag_int.txt");
	//file_drag_int << drag_int;
	//file_drag_int.close();

	//int period = T;
	//ofstream file_PF;
	//file_PF.open("../examples/period.txt");
	//file_PF << period;
	//file_PF.close();
	//cout << "period: " << period << endl;

	int PhaseCal_flag = 0;
	// save PhaseCal_flag
	ofstream file_PF_flag;
	file_PF_flag.open("../examples/PhaseCal_flag.txt");
	file_PF_flag << PhaseCal_flag;
	file_PF_flag.close();
	cout << "PhaseCal_flag: " << PhaseCal_flag << endl;
	//

	// Set plate motion
	plate.setMotion(motion);
	Geometry geom;
	geom.addBody(plate);
	//geom.moveBodies(0);

	// Setup equations to solve

	double alpha = 0;  // angle of background flow
	BaseFlow q_potential(grid, magnitude, alpha);

	// Transfer the motion from inertial coordinate to body coordinate
	Motion* m = geom.transferMotion(); // pull motion from first RigidBody object
	geom.transferCenter(xC, yC);  // pull center of motion from RigidBody object
	q_potential.setMotion(*m);
	q_potential.setCenter(xC, yC);
	// Setup IBPM Solver
	cout << "Setting up Navier Stokes model..." << flush;
	NavierStokesModel model(grid, geom, Reynolds, q_potential);
	model.init();
	cout << "done" << endl;
	//cout << "amplitude = " << amplitude << endl;

	// Setup timestepper
	// NonlinearIBSolver solver(grid, model, dt, Scheme::EULER);
	// NonlinearIBSolver solver(grid, model, dt, Scheme::AB2);
	NonlinearIBSolver solver(grid, model, dt, Scheme::RK3);
	solver.init();
	//
	State x(grid, geom.getNumPoints());
	x.omega = 0.1;
	x.q = 0.1;
	x.f = 0.1;

	// if not the first ieration, load bin file that gives better initial gauess
	if (k_global != 0) {
		x.load(outdir_PBF + "pitch000.bin");
	}

	State Boundary(grid, geom.getNumPoints());
	Boundary.f = 0.1;
	//


	//double cycle = 1 / freq / dt;
	//int cycle_int = (int)cycle;
//////////////////////////////
	int T = size1;
	int numSteps = (int)T * Ini_step;

	///////////////////////////////
//////////////////////////////
///////////////////////////////
// for test only
	double drag_save[numSteps];
	double lift_save[numSteps];
	////////////////////////////////
	////////////////////////////////
	////////////////////////////////

	for (int i = 1; i <= numSteps; ++i) {
		// double theta = 1* sin(2 * PI / T * x.time);


		solver.advance(x);
		double xF, yF;
		x.computeNetForce(xF, yF);
		double theta = motion.gettheta(x.time);
		double thetadot = motion.getthetadot(x.time);
		double ydot = motion.getydot(x.time);
		double y = motion.gety(x.time);

		q_potential.setAlphaMag(x.time);
		alpha = q_potential.getAlpha();
		drag = xF * cos(theta) - yF * sin(theta);
		lift = xF * sin(theta) + yF * cos(theta);

		////////
		cout << " step " << setw(4) << i << " time " << setw(5) << x.time << endl;;
		cout << " lift ..." << lift << " drag ..." << drag << endl;
		// cout << " theta ... " << theta << endl;
		// cout << " thetadot ... " << thetadot << endl;
		// cout << " ydot ... " << ydot << endl;
		// cout << " y ... " << y << endl;
		/////////////////////////////
		Boundary.f = plate.getPoints();

		double zM = 0.0;
		double drag_check = 0.0;
		double lift_check = 0.0;
		for (int j = 0;j < geom.getNumPoints();++j) {
			//zM = x.f(0, j) * grid.Dx() * grid.Dx() * (Boundary.f(1, j) - yC) * sin(theta) + x.f(1, j) * grid.Dx() * grid.Dx() * (Boundary.f(0, j) - xC) * cos(theta) + zM;
			zM = x.f(1, j) * grid.Dx() * grid.Dx() * (Boundary.f(0, j) - xC) + zM;
			drag_check = x.f(0, j) * grid.Dx() * grid.Dx() + drag_check;
			lift_check = x.f(1, j) * grid.Dx() * grid.Dx() + lift_check;
		}

		drag_save[i - 1] = drag;
		lift_save[i - 1] = lift;
	}
	///////////////////////////
	//////////////////////////
	// for test only
	// Save dDdt
	ofstream myfile_eta("../examples/" + outdir + "D_save.txt");
	if (myfile_eta.is_open())
	{
		// for (int count = 0; count < size; count++) {
		for (int count = 0; count < numSteps; ++count) {
			myfile_eta << drag_save[count] << " ";
		}
		myfile_eta.close();
	}
	else cout << "Unable to open file";

	// Save dDdt
	ofstream myfile_zeta("../examples/" + outdir + "L_save.txt");
	if (myfile_zeta.is_open())
	{
		// for (int count = 0; count < size; count++) {
		for (int count = 0; count < numSteps; ++count) {
			myfile_zeta << lift_save[count] << " ";
		}
		myfile_eta.close();
	}
	else cout << "Unable to open file";
	////////////////////////
	/////////////////////////


	// quiet nan

	//string icfile = "pitch200.bin";
	//x.load(icfile);
	//cout << "Loading initial condition from file: " << icfile << endl;
	//

	// Convert Vec variables to State
	int M = 0;
	//Vec Omega_in(((nx-1) * (ny-1) * 3) * ngrid + geom.getNumPoints() * 2);
	Vec Omega_in(((nx - 1)* (ny - 1) + (nx + 1) * ny + nx * (ny + 1))* ngrid + geom.getNumPoints() * 2);

	for (int k = 0; k < ngrid; ++k){
		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				Omega_in(M) = x.omega(k, i, j);
				M = M + 1;
			}

		}

		for (int i = 0; i <= nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				Omega_in(M) = x.q(k, 0, i, j);
				M = M + 1;
			}

		}

		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j <= ny; ++j) {
				Omega_in(M) = x.q(k, 1, i, j);
				M = M + 1;
			}

		}

	}


	for (int i = 0; i < geom.getNumPoints(); ++i) {
		Omega_in(M) = x.f(0, i);
		M = M + 1;
	}

	for (int i = 0; i < geom.getNumPoints(); ++i) {
		Omega_in(M) = x.f(1, i);
		M = M + 1;
	}


	//for (int i = 0; i < _numPoints; ++i) {
	//	Omega(M) = Fx.f(0, i);
	//	M = M + 1;
	//}

	//for (int i = 0; i < _numPoints; ++i) {
	//	Omega(M) = Fx.f(1, i);
	//	M = M + 1;
	//}
	// Newton-GMRES solver
	linear_solver_ns::Gmres<Vec> gmres(innerproduct, 5, 1.e-3); // 1.e-3
	// constexpr double tol{ Newton_tol };
	constexpr double tol{ 1.e-3 };  // 1.e-8  1.e-3
	constexpr int max_iter{ 15 };
	constexpr double jacobian_dx{ 1.e-3 }; // 1.e-3
	std::unique_ptr<newton_ns::Newton<Vec>> newton;
		newton.reset(
			new newton_ns::Newton<Vec>{ F, gmres, norm, tol, max_iter,
					jacobian_dx, 1, true }
		);

	//newton.reset(
	//	new newton_ns::NewtonArmijo<Vec>{ F, gmres, norm, tol, max_iter,
	//			jacobian_dx, 1, 10, 0.1, 0.5, 1.e-4, true }
	//);

	try {
		newton->solve(Omega_in);
	}

	catch (newton_ns::NewtonError& ne) {
		cerr << "Newton error: " << ne.what() << endl;
		exit(2);
	}
	catch (std::exception& e) {
		cerr << "Other exception: " << e.what() << endl;
		exit(3);
	}

	// Load drag intergral
	//std::ifstream input_file_drag_int("../examples/drag_int.txt");
	////double drag_int;
	//input_file_drag_int >> drag_int;
//////////////
  save_flag = 1;
  F(Omega_in);
	save_flag = 0;
/////////////
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	// End timing
	std::cout << "It took me " << time_span.count() << " seconds.";
	std::cout << std::endl;

	if (cost == 0){
		drag_int = CT;
	}
	else{
		drag_int = CT/CP;
	}


	cout << "drag_int: " << drag_int << endl;
	return drag_int;
}

/////////////
// inline double drag_int_fun(const Vec& x) {
// 	//cout << "x0..." << x(0) << endl;
// 	//cout << "x1..." << x(1) << endl;
// 	 double y = ((x(0) - 3) * (x(0) - 3)) + ((x(1) + 5) * (x(1) + 5));
// 	 //cout << "y..." << y << endl;
// 	// double y = (x(0) - 3) * (x(0) - 3);
// 	return y;
// }
// compute periodic orbit for ADJ Solver

inline Vec ADJ_F(const Vec& ADJ_x_in){
	double lift = 0.;
	double drag = 0.;
	Grid grid( nx, ny, ngrid, length, xOffset, yOffset);
	// Make a flat plate, length 1, with center at 1/4 chord
	RigidBody plate;
	plate.addLine(0, 0, 1, 0, grid.Dx());
	plate.setCenter(0, 0);

	//
	int TT = size1;
	int end_index = TT * num_cycles;  // only use end_index = TT for test
	// double pitchPeriod_temp = round(TT / dt);
	// // TT = pitchPeriod_temp * dt;
	// int pitchPeriod = (int)pitchPeriod_temp;

	double pitch_in[size1*num_cycles+1];
	double plunge_in[size1*num_cycles+1];
	double motiontime_in[size1*num_cycles+1];


	for(int i = 0; i<size1*num_cycles+1; ++i){
				// pitch_in[i] =  input(i,0);
				// plunge_in[i] =  input(i,1);
				// motiontime_in[i] = MotionTime(i);

				pitch_in[i] =  MotionParameter(i,0);
				plunge_in[i] =  MotionParameter(i,1);
				motiontime_in[i] = MotionTime(i);


	}

	// create dqdb
	std::ifstream file_tc("../examples/theta_adj.txt");
	double t_t[size1*num_cycles+1];
	// Update Error
	for (int i = 0; i < size1*num_cycles+1; ++i)
	 {
		 file_tc >> t_t[i];
	 }
	// Save dqdb
	ofstream myfile_thetacheck("../examples/theta_adj.txt");
	if (myfile_thetacheck.is_open())
	{
	 // for (int count = 0; count < size; count++) {
	 for (int count = 0; count <size1*num_cycles+1; ++count) {
		 myfile_thetacheck << pitch_in[count] << " ";
	 }
	 myfile_thetacheck.close();
	}
	else cout << "Unable to open file";
	////////////

	PitchPlungeT_mod motion(
		dt,
		motiontime_in,
		pitch_in,
		plunge_in,
		end_index
		// matrix.row(i)
		// matrix.col(j)

	);
	plate.setMotion(motion);
	Geometry geom;
	geom.addBody(plate);
	// geom.moveBodies(0);

	// Setup equations to solve

	double alpha = 0;  // angle of background flow
	int periodStart = 0;
	int period = TT;
	BaseFlow q_potential(grid, magnitude, alpha);

/////////////////////////////
	// Motion* mot = geom.transferMotion(); // pull motion from first RigidBody object
	// geom.transferCenter(xC, yC);  // pull center of motion from RigidBody object
	// q_potential.setMotion(*mot);
	// q_potential.setCenter(xC, yC);
//////////////////////////
	NavierStokesModel* model = NULL;
	AdjointIBSolver2* adjointsolver = NULL;
	IBSolver* solver = NULL;

	/*State x00(grid, geom.getNumPoints());
	vector<State> x0(period, x00);*/
/////////////////
	Motion* mot = NULL; // this is for adjoint tests; get rid of later
///////////////
/////////////////////////////
	mot = geom.transferMotion(); // pull motion from first RigidBody object
	geom.transferCenter(xC, yC);  // pull center of motion from RigidBody object
	q_potential.setMotion(*mot);
	q_potential.setCenter(xC, yC);
//////////////////////////

	//for (int i = 0; i < period; i++){  // if change period here, also change the period in vector<State> x0(period, x00);
	//	// sprintf("PBF", "PBF".c_str(), i + periodStart);
	//	// if ( ! x0[i].load(pbffilename) ) {

	//	cout << "i = " << std::to_string(i) << endl;

	//	if (i < 10) {
	//		cout << "first" << endl;
	//		if ( ! x0[i].load(outdir_PBF + "pitch00" + std::to_string(i) + ".bin") ) {
	//				cout << "base flow " << "PBF/pitch00" + std::to_string(i) + ".bin" << " failed to load. Exiting program." << endl;
	//				exit(1);
	//			 }

	//	}
	//	if(i >= 10 && i < 100){
	//		cout << "second" << endl;
	//		if ( ! x0[i].load(outdir_PBF + "pitch0" + std::to_string(i) + ".bin") ) {
	//				cout << "base flow " << "PBF/pitch0" + std::to_string(i) + ".bin" << " failed to load. Exiting program." << endl;
	//				exit(1);
	//			 }

	//	}
	//	if(i >= 100){
	//		cout << "third" << endl;
	//		if ( ! x0[i].load(outdir_PBF + "pitch" + std::to_string(i) + ".bin") ) {
	//				cout << "base flow " << "PBF/pitch" + std::to_string(i) + ".bin" << " failed to load. Exiting program." << endl;
	//				exit(1);
	//			 }

	//	}

	//}


	//x00 = x0[0];
	///////////////////////
	// mot = geom.transferMotion();
	////////////////////
	model = new NavierStokesModel( grid, geom, Reynolds);
	// model = new NavierStokesModel( grid, geom, Reynolds, q_potential );


	Scheme::SchemeType str2scheme( string integratorType );
	string integratorType = "rk3";
	// string integratorType = "AB2"; //////////////////////
	// string integratorType = "EULER"; //////////////////////
	Scheme::SchemeType schemeType = str2scheme( integratorType );

	adjointsolver = new AdjointIBSolver2( grid, *model, dt, schemeType, x0_glob, period, cost, CT, CP, *mot );
	// adjointsolver = new AdjointIBSolver2( grid, *model, dt, "rk3", x0, period, *mot );
	solver = adjointsolver;

	State x( grid, geom.getNumPoints() );
	x.omega = 0.;
	x.f = 0.;
	x.q = 0.;

	State Boundary(grid, geom.getNumPoints());
	Boundary.f = 0.1;
		// update the geometry to the current time
	// geom.moveBodies( x.time );

	// Initialize model and timestepper
	model->init();
	cout << "using " << solver->getName() << " timestepper" << endl;
	cout << "    dt = " << dt << "\n" << endl;

			solver->init();
			solver->save( "outdir_name" );

	// Calculate flux for state, in case only vorticity was saved
	if( ! q_potential.isStationary() ) {
			q_potential.setAlphaMag( x.time );
			alpha = q_potential.getAlpha();
			cout << "q_potential calculation" << endl;
	}
	// model->updateOperators( x.time ); /////////////////////////
	// model->refreshState( x );  ////////////////////////

	cout << endl << "Initial timestep = " << x.timestep << "\n" << endl;


	StateVector xv(x);
////////////
// pass the initial value to the state
	int N = 0;

	for (int k = 0; k<ngrid;++k){
		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				xv.x.omega(k, i, j) = ADJ_x_in(N);   //x.omega(0,i,j)
				N = N + 1;
			}

		}

		for (int i = 0; i <= nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				xv.x.q(k, 0, i, j) = ADJ_x_in(N);   //x.omega(0,i,j)
				N = N + 1;
			}

		}

		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j <= ny; ++j) {
				xv.x.q(k, 1, i, j) = ADJ_x_in(N);   //x.omega(0,i,j)
				N = N + 1;
			}

		}

	}


	for (int i = 0; i < geom.getNumPoints(); ++i) {
		xv.x.f(0, i) = ADJ_x_in(N);   //x.omega(0,i,j)
		N = N + 1;
	}

	for (int i = 0; i < geom.getNumPoints(); ++i) {
		xv.x.f(1, i) = ADJ_x_in(N);   //x.omega(0,i,j)
		N = N + 1;
	}

 ////////////////////
	int numSteps = num_cycles * period;

	 // start adjoint solver
	 Vec objgrad;
	 int k;
	 objgrad.resize(numSteps+1); // vector to store gradient
	 double dx2 = grid.Dx() * grid.Dx();
	 double dx = grid.Dx();
	 double alpha_previous;
	 double beta_previous;
	 double alpha_adj;
	 double beta_adj;


	 Vec dA_check(size1);
	 Vec increment(size1);
	 Vec increment1(size1);

	 double lift_adj = 0.;
	 double drag_adj = 0.;
	 double drag_adj_pre = 0.;
	 double y = 0.;
	 double y_pre = 0.;
   double dgdT;

	 Vec dgdtheta_save(numSteps+1);
	 Vec dDdy(numSteps+1);
	 Vec dDdtheta(numSteps+1);
	 Vec eta(numSteps+1);
	 Vec zeta(numSteps+1);

	 double dqdaN_s[numSteps+1];

	 double dqdbN_s[numSteps+1];

	 double dudaphi_s[numSteps+1];

	 double dudbphi_s[numSteps+1];

	 double theta_adj_check[numSteps+1];

	 // Create output directory, if does not already exist
	mkdir(outdir_ADJ.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

	// Setup output routines
	// if (save_flag == 1){
	OutputForce force(outdir_ADJ + "ADJ_force.dat");
	OutputTecplot tecplot(outdir_ADJ + "ADJ%03d.plt", "Pitching plate, step %03d", 1);
	if(TecplotAllGrids) tecplot.setFilename(outdir_ADJ + "ADJ%03d_g%01d.plt" );

	// }

	OutputRestart restart(outdir_ADJ + "ADJ%03d.bin");
	Logger logger;
	//Output Tecplot file every few timesteps

	if (save_flag == 1){
		numSteps = num_cycles * period + 1;
		xv.x.load(outdir_ADJ + "ADJ000.bin");
		logger.addOutput(&tecplot, 10 );
		logger.addOutput(&force, 10 );
		logger.addOutput(&restart, 1000);
	}
	else{
		numSteps = period;
		logger.addOutput(&restart, 1000);
	}


	logger.init();
//	logger.doOutput(xv.x);


	double zm_ini = 0;
	double theta_ini = motion.gettheta(0 * dt);

	// logger.doOutput(theta_ini, q_potential, xv.x, zm_ini);
  logger.doOutput(xv.x);
//////////////////////
// tranform the flow from NS using unsteady base flow to steady base flow !!!!
//////////////////////
	// for(int i = 0; i <= numSteps; i++) {
	for(int i = 0; i < numSteps; i++) {
	 // for(int i = 1; i <= numSteps; i++) {

		 // if(i > period){
			//  k = std::remainder(i,period);
		 // }
		 // if(k < 0){
			//  k = k + period;
		 // }


		 	 //solver->advance(xv.x);  ///// maybe put it at the beginning of the loop
			 //k = period - xv.x.timestep - 1; // if not periodic
			 k = (period - (xv.x.timestep % period)) % period; // if periodic  have advance(x) before this line?
			 double partial = 0.;
			 Flux dqdp(x0_glob[k].q);
			 TangentSE2 gg(0, 0, 0, 0, 1., 0);  // eq 14 in Daniel's note, add the one for eq 15 which is the second argument.
			 // (x, y, theta, xdot, ydot, thetadot).
			 dqdp.setFlow(gg, xC, yC);
			 double deriv1 = (-InnerProduct(dqdp, CrossProduct(xv.x.q, x0_glob[k].omega))) / dx2;
			 BoundaryVector dudp(xv.x.f.getNumPoints());
			 dudp = -model->flux2boundary(dqdp);
			 double deriv2 = InnerProduct(dudp, xv.x.f);

			 // logger.doOutput( q_potential, x);

			 objgrad(i) = partial + deriv1 + deriv2;


			 y_pre = y;
			 y = motion.gety(k * dt);
			 double theta = motion.gettheta(k * dt);
			 double thetadot = motion.getthetadot(k * dt);
			 double ydot = motion.getydot(k * dt);
			 double y = motion.gety(k * dt);


			 cout << " k ... " << k << endl;
			 cout << " x.timestep ... " << xv.x.timestep << " x.time ... " << xv.x.time << endl;
			 // cout << " theta_pass ... " << theta << endl;
			 // cout << " thetadot_pass ... " << thetadot << endl;
			 // cout << " ydot_pass ... " << ydot << endl;
			 cout << "line... 1429" << endl;

			 ///// dgd()
			 double dgdh = 0.;
			 double dgdhdot = 0.;
			 double dgdtheta = 0.;
			 double dgdthetadot = 0.;
			 dgdT = 0.;
			 double xF, yF, zM;
			 x0_glob[k].computeNetForce(xF, yF);
			 // to check if the time step is right or note can be removed
			 drag_adj_pre = drag_adj;
			 drag_adj = xF * cos(theta) - yF * sin(theta);
			 lift_adj = xF * sin(theta) + yF * cos(theta);
			 cout << "line... 1443" << endl;
			 // compute the pitching moment zM
			 Boundary.f = plate.getPoints();
             zM = 0.0;
			 for (int j = 0;j < geom.getNumPoints();++j) {
			 //zM = x.f(0, j) * grid.Dx() * grid.Dx() * (Boundary.f(1, j) - yC) * sin(theta) + x.f(1, j) * grid.Dx() * grid.Dx() * (Boundary.f(0, j) - xC) * cos(theta) + zM;
			 // zM = x.f(1, j) * grid.Dx() * grid.Dx() * (Boundary.f(0, j) - xC) + zM;
			 zM = x0_glob[k].f(1, j) * grid.Dx() * grid.Dx() * (Boundary.f(0, j) - xC) + zM;
			 // cout << " x0[k].f(1, j) ... " << x0[k].f(1, j) << endl;
			 // cout << " grid.Dx() ... " << grid.Dx() << endl;
			 // cout << " (Boundary.f(0, j) - xC) ... " << (Boundary.f(0, j) - xC) << endl;
			 // cout << " zM ... " << zM << endl;
			 }
			 cout << "line... 1456" << endl;
			 if (cost == 0){
				 // cost = 0, thrust
			 	// dgdtheta = -(-xF * sin(theta) - yF * cos(theta))*2; // * 2
			 	dgdtheta = (-xF * sin(theta) - yF * cos(theta)) * 2 / ((double)period*dt) * magnitude; // * 2  add /(period * dt) for everything related to g
			 	// dgdtheta = (xF * sin(-theta) + yF * cos(-theta))*2; // * 2
                dgdtheta_save(i) = dgdtheta;

			 	dgdthetadot = 0 / ((double)period*dt);
			 	dgdh = 0 / ((double)period*dt);
			 	dgdhdot = 0 / ((double)period*dt);
			 	// dgdT = (-1) * pow((double)period*dt,-2.) * CT * 2;
				dgdT = (-1) * CT / ((double)period * dt);
				// cout << "(double)period ... " << (double)period << endl;
				// cout << "(double)period*dt ... " << (double)period*dt << endl;
				// cout << "pow((double)period*dt,-2) ... " << pow((double)period*dt,-2) << endl;
				// cout << "CT ... " << CT << endl;
			 }
			 else if (cost == 1){
				 cout << "line... 1474" << endl;
				 // cost = 1, efficiency
				 dgdtheta = 1/CP * (- 2 * sin(theta) *  xF - 2 * cos(theta) * yF) / ((double)period * dt) - CT/(pow(CP,2)) * (2 * sin(theta) * yF - 2 * cos(theta) * xF) * ydot / ((double)period * dt);
         // dgdtheta = (2 * sin(theta) * yF * ydot - 2 * cos(theta) * xF * ydot - 2 * sin(theta) * xF - 2 * cos(theta) * yF) / ((double)period*dt);
				 // dgdtheta = (-2 * sin(theta) * yF * ydot + 2 * cos(theta) * xF * ydot - 2 * sin(theta) * xF - 2 * cos(theta) * yF);

				 cout << "line... 1481" << endl;
				 // dgdthetadot = -(-zM);  // -
				 dgdthetadot = -(-(CT/(pow(CP,2)) * (-2 * zM) / ((double)period * dt)));
				 // dgdthetadot = -(-2 * zM) / ((double)period*dt);  // -
                cout << " zM ... " << zM << endl;

				 dgdh = 0 / ((double)period*dt);

				 dgdhdot = -(-(CT/(pow(CP,2)) * (-2 *cos(theta) * yF - 2 * sin(theta) * xF) / ((double)period * dt)));
				 // dgdhdot = -(-(CT/(pow(CP,2)) * (-2 *cos(theta) * yF - sin(theta) * xF)));
				 // dgdhdot = -(-2 * cos(theta) * yF -2 * sin(theta) * xF) / ((double)period*dt);  /// -

				 dgdT = 0;
				 // dgdT = (-1) * pow((double)period*dt,-2) * (CP + CT) * dt * 2;
				 cout << "line... 1495" << endl;

			 	// // Note that the dgd() in the note is -dgd()
			 	// dgdtheta = (-1) * ((xF * sin(theta) + yF * cos(theta)) * (-cos(theta) * yF * ydot -sin(theta) * xF * ydot - zM * thetadot)/(pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period\
			 	// 					- (-cos(theta) * xF + yF * sin(theta)) * (sin(theta) * yF * ydot -cos(theta) * xF * ydot)/(pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period);
				//
			 	// dgdthetadot = (-1) * (0/(pow(-cos(theta) * yF *ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period\
			 	// 						 - (-xF * cos(theta) + yF * sin(theta) * (-1) * zM)/(pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period);
				//
			 	// dgdh = 0;
				//
			 	// dgdhdot	= (-1) * (0/(pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period\
			 	// 					 - (-xF * cos(theta) + yF * sin(theta)) * (-cos(theta) * yF - sin(theta) * xF) / (pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period);
				//
			 	// dgdT = (-1) * pow(period,-2) * xF/(yF * ydot + zM * thetadot);
			 }
			 //////////////////

// alpha
State x_dqda(grid, geom.getNumPoints());
Flux dqda(x_dqda.q);

Flux::index ind;
for (int lev=0; lev<dqda.Ngrid(); ++lev) {

		double xdot_dqda = sin(theta); //-
		double ydot_dqda = cos(theta); //-
		// double xdot = 0;
		// double ydot = 1;
		double xdiff, ydiff;


		// double dx = dqda.Dx(lev);
		for(ind = dqda.begin(X); ind != dqda.end(X); ++ind) {
				xdiff = dqda.x(lev,ind) - xC;
				ydiff = dqda.y(lev,ind) - yC;
				// dqda(lev,ind) = (xdot -thetadot*ydiff)*dx;
				dqda(lev,ind) = (xdot_dqda)*dx;
		}
		for(ind = dqda.begin(Y); ind != dqda.end(Y); ++ind) {
				xdiff = dqda.x(lev,ind) - xC;
				ydiff = dqda.y(lev,ind) - yC;
				// dqda(lev,ind) = (ydot + thetadot*xdiff)*dx;
				dqda(lev,ind) = (ydot_dqda)*dx;
		}
}
cout << "line... 1541" << endl;

// double deriv1_alpha = (-InnerProduct(dqda, CrossProduct(xv.x.q, x0[k].omega))) / dx2;  // cout this term and the one in next line
double deriv1_alpha = (InnerProduct(dqda, CrossProduct(xv.x.q, x0_glob[k].omega))) / dx2;  // cout this term and the one in next line ///////mod
BoundaryVector duda(xv.x.f.getNumPoints());
duda = -model->flux2boundary(dqda); // remove -
// xv.x.f = -1 * xv.x.f;   ////////////////////////////
double deriv2_alpha = -InnerProduct(duda, xv.x.f);
alpha_previous = alpha_adj;
alpha_adj = dgdhdot + deriv1_alpha + deriv2_alpha;
double dalphadt = (alpha_adj - alpha_previous) / dt;  /// add -  !!!

cout << "line... 1553" << endl;
////////////////////////////////
////////
// cout << " deriv1_alpha ... " << deriv1_alpha << endl;
// cout << "deriv2_alpha ... " << deriv2_alpha << endl;
/////////////////
State x_dqdb(grid, geom.getNumPoints());
Flux dqdb(x_dqdb.q);
Flux::index ind_1;
for (int lev=0; lev<dqdb.Ngrid(); ++lev) {

		double xdiff, ydiff;

		// double dx = dqda.Dx(lev);
		for(ind_1 = dqdb.begin(X); ind_1 != dqdb.end(X); ++ind_1) {
				xdiff = dqdb.x(lev,ind_1) - xC;
				ydiff = dqdb.y(lev,ind_1) - yC;
				// dqda(lev,ind) = (xdot -thetadot*ydiff)*dx;
				dqdb(lev,ind_1) = -(ydiff)*dx;
		}
		for(ind_1 = dqdb.begin(Y); ind_1 != dqdb.end(Y); ++ind_1) {
				xdiff = dqdb.x(lev,ind_1) - xC;
				ydiff = dqdb.y(lev,ind_1) - yC;
				// dqda(lev,ind) = (ydot + thetadot*xdiff)*dx;
				dqdb(lev,ind_1) = -(-xdiff)*dx;
		}
}

			 // double deriv1_beta = (-InnerProduct(dqdb, CrossProduct(xv.x.q, x0[k].omega))) / dx2;
			 double deriv1_beta = (InnerProduct(dqdb, CrossProduct(xv.x.q, x0_glob[k].omega))) / dx2;  //////////////////////// mod
			 BoundaryVector dudb(xv.x.f.getNumPoints());
			 dudb = -model->flux2boundary(dqdb);
			 double deriv2_beta = -InnerProduct(dudb, xv.x.f);
			 beta_previous = beta_adj;

			 beta_adj = dgdthetadot + deriv1_beta + deriv2_beta;
			 double dbetadt = (beta_adj - beta_previous) / dt;
			 /////////////////
			 // cout << " dalpha/dt ..." << dalphadt << " dbeta/dt ..." << dbetadt << endl;
			 // cout << " dgdthetadot ..." << dgdthetadot << " beta1 ..." << deriv1_beta << " beta2 ..." << deriv2_beta << endl;

///////////
State x_dqdh(grid, geom.getNumPoints());
Flux dqdh(x_dqdb.q);
Flux::index ind_2;
for (int lev=0; lev<dqdh.Ngrid(); ++lev) {

		double xdiff, ydiff;

		// double dx = dqda.Dx(lev);
		for(ind_2 = dqdh.begin(X); ind_2 != dqdh.end(X); ++ind_2) {
				xdiff = dqdh.x(lev,ind_2) - xC;
				ydiff = dqdh.y(lev,ind_2) - yC;

				dqdh(lev,ind_2) = -thetadot * sin(theta) * dx * 0;  //////////////////////
		}
		for(ind_2 = dqdh.begin(Y); ind_2 != dqdh.end(Y); ++ind_2) {
				xdiff = dqdh.x(lev,ind_2) - xC;
				ydiff = dqdh.y(lev,ind_2) - yC;

				dqdh(lev,ind_2) = thetadot * cos(theta) * dx * 0;  //////////////////
		}
}

BoundaryVector dudh(xv.x.f.getNumPoints());
dudh = -model->flux2boundary(dqdh);
///////////
			 // Flux dqdh(x0[k].q);
			 // dqdh = dqdh * 0;
			 // BoundaryVector dudh(xv.x.f.getNumPoints());
			 // dudh = -model->flux2boundary(dqdh);
////////////





/////////////////////
/////////////////
State x_dqdtheta(grid, geom.getNumPoints());
Flux dqdtheta(x_dqdtheta.q);
Flux::index ind_3;
for (int lev=0; lev<dqdtheta.Ngrid(); ++lev) {

		double xdiff, ydiff;

		// double dx = dqda.Dx(lev);
		for(ind_3 = dqdtheta.begin(X); ind_3 != dqdtheta.end(X); ++ind_3) {
				xdiff = dqdtheta.x(lev,ind_3) - xC;
				ydiff = dqdtheta.y(lev,ind_3) - yC;
				// dqdtheta(lev,ind_3) = (-sin(theta) - thetadot*xdiff  -ydot*cos(theta))*dx;
				dqdtheta(lev,ind_3) = (-sin(theta)-ydot*cos(theta))*dx;  ////////////////////////
		}
		for(ind_3 = dqdtheta.begin(Y); ind_3 != dqdtheta.end(Y); ++ind_3) {
				xdiff = dqdtheta.x(lev,ind_3) - xC;
				ydiff = dqdtheta.y(lev,ind_3) - yC;
				// dqdtheta(lev,ind_3) = (-cos(theta) + thetadot*ydiff +ydot*sin(theta))*dx;
				dqdtheta(lev,ind_3) = (-cos(theta)+ydot*sin(theta))*dx;   /////////////////////////
		}
}
/////////////////////
			 // Flux dqdtheta(x0[k].q);
			 // TangentSE2 gg_theta(0, 0, 0, -sin(theta)-ydot*cos(theta) , -cos(theta)+ydot*sin(theta), 0);  // eq 14 in Daniel's note, add the one for eq 15 which is the second argument. The right one
			 // // (x, y, theta, xdot, ydot, thetadot).
			 // dqdtheta.setFlow(gg_theta, xC, yC);
///////////////////
			 BoundaryVector dudtheta(xv.x.f.getNumPoints());
			 dudtheta = -model->flux2boundary(dqdtheta);

			 double delta_h;
			 double delta_theta;
			 // delta_h = dgdh-(InnerProduct(dqdh, CrossProduct(xv.x.q, x0[k].omega)))/dx2-InnerProduct(dudh, xv.x.f)-dalphadt; //
			 // x0[k].omega = x0[k].omega + 2 * thetadot;///////////////////////////////////////////////////////
			 delta_h = dgdh+(InnerProduct(dqdh, CrossProduct(xv.x.q, x0_glob[k].omega)))/dx2-InnerProduct(dudh, xv.x.f)-dalphadt; // ////////////// mod
			 // delta_theta = dgdtheta-(InnerProduct(dqdtheta, CrossProduct( xv.x.q, x0[k].omega)))/dx2-InnerProduct(dudtheta, xv.x.f)-dbetadt;
			 delta_theta = dgdtheta+(InnerProduct(dqdtheta, CrossProduct( xv.x.q, x0_glob[k].omega)))/dx2-InnerProduct(dudtheta, xv.x.f)-dbetadt; ////////////// mod
			 cout << "line... 1669" << endl;
			 cout << " delta_h ..." << delta_h << " delta_theta ..." << delta_theta << endl;
			 //////////
			 // cout << " dgdh ... " << dgdh << endl;
			 // double dqdaN = (InnerProduct(dqdh, CrossProduct(xv.x.q, x0[k].omega)))/dx2;
			 // cout << " dqda*N ... " << dqdaN << endl;
			 //
			 // double dqdbN = (InnerProduct(dqdtheta, CrossProduct(xv.x.q, x0[k].omega)))/dx2;
			 // cout << " dqdb*N ... " << dqdbN << endl;
			 //
			 // double dudaphi = InnerProduct(dudh, xv.x.f);
			 // cout << " duda*phi ... " << dudaphi << endl;
			 // double dudbphi = InnerProduct(dudtheta, xv.x.f);
			 // cout << " dudb*phi ... " << dudbphi << endl;
			 // cout << " dalphadt ... " << dalphadt << endl;
			 // double Ssin = sin(xv.x.time/2*(2*PI));
			 // cout << " sin(t) ... " << Ssin << endl;

			 //////////
			 eta(i) = delta_h;
			 zeta(i) = delta_theta;
      // dA amp for heaving
			// dB phase for heaving
			// dC amp for theta
			// dD phase for theta

			// for (int n = 0; n < size3; ++n){
			//   /// add dgdA ...
			// 		dA(n) = dA(n) + eta(i) * sin(2 * PI * (n+1) * xv.x.time/((double)period*dt) + B(n)) * dt;   /// keep dt or not
			// 		dB(n) = dB(n) + eta(i) * A(n) * cos(2 * PI * (n+1) * xv.x.time/((double)period*dt) + B(n)) * dt;
			// 		dC(n) = dC(n) + zeta(i) * sin(2 * PI * (n+1) * xv.x.time/((double)period*dt) + D(n)) * dt;
			// 		dD(n) = dD(n) + zeta(i) * C(n) * cos(2 * PI * (n+1) * xv.x.time/((double)period*dt) + D(n)) * dt;
			//
			// }

			// eta(i) = delta_h;
			// zeta(i) = delta_theta;
			// for (int n = 0; n < size3; ++n){
			// 	/// add dgdA ...
			// 		dA(n) = dA(n) + eta(i) * sin(2 * PI * (n+1) * (double)k*dt/((double)period*dt) + B(n)) * dt;   /// keep dt or not
			// 		dB(n) = dB(n) + eta(i) * A(n) * cos(2 * PI * (n+1) * (double)k*dt/((double)period*dt) + B(n)) * dt;
			// 		dC(n) = dC(n) + zeta(i) * sin(2 * PI * (n+1) * (double)k*dt/((double)period*dt) + D(n)) * dt;
			// 		dD(n) = dD(n) + zeta(i) * C(n) * cos(2 * PI * (n+1) * (double)k*dt/((double)period*dt) + D(n)) * dt;
			//
			// 		dA_check(i) = dA(size3-1);
			// 		increment(i) = eta(i) * sin(2 * PI * (size3) * (double)k*dt/((double)period*dt) + B(size3-1)) * dt;
			// 		increment1(i) = eta(i) * sin(2 * PI * (size3) * (double)k*dt/((double)period*dt) + B(size3-1));
			//
			// }


			 // if (i == 0){
				//  dDdy(i) = delta_h;
				//  dDdtheta(i) = delta_theta;
			 // }
			 // else{
				//  dDdy(i) = eta * (sin(2*PI));
				//  dDdtheta(i) = delta_theta;
			 // }




			 // dqdaN_s[i] = dqdaN;
			 // dqdbN_s[i] = dqdbN;
			 // dudaphi_s[i] = dudaphi;
			 // dudbphi_s[i] = dudbphi;
       theta_adj_check[i] = theta;

			 cout << " dDdy ..." << dDdy(i) << endl;

       double zm = 0;
	   solver->advance(xv.x);  ///// maybe put it at the beginning of the loop
       logger.doOutput(xv.x);
			 // logger.doOutput( theta, q_potential, xv.x, zm);

	 }


	 //////////////////////////////
	 // Convert State variables in to Vec
	 int M = 0;
	 //Vec Omega(((nx-1) * (ny-1) * 3) * ngrid + geom.getNumPoints() * 2);
	 Vec Omega(((nx - 1)* (ny - 1) + (nx + 1) * ny + nx * (ny + 1))* ngrid + geom.getNumPoints() * 2);

	 for (int k = 0; k < ngrid; ++k){
	 	for (int i = 1; i < nx; ++i) {
	 		for (int j = 1; j < ny; ++j) {
	 			Omega(M) = xv.x.omega(k, i, j);
	 			M = M + 1;
	 		}

	 	}

	 	for (int i = 0; i <= nx; ++i) {
	 		for (int j = 0; j < ny; ++j) {
	 			Omega(M) = xv.x.q(k, 0, i, j);
	 			M = M + 1;
	 		}

	 	}

	 	for (int i = 0; i < nx; ++i) {
	 		for (int j = 0; j <= ny; ++j) {
	 			Omega(M) = xv.x.q(k, 1, i, j);
	 			M = M + 1;
	 		}

	 	}

	 }


	 for (int i = 0; i < geom.getNumPoints(); ++i) {
	 	Omega(M) = xv.x.f(0, i);
	 	M = M + 1;
	 }

	 for (int i = 0; i < geom.getNumPoints(); ++i) {
	 	Omega(M) = xv.x.f(1, i);
	 	M = M + 1;
	 }

	 double f_norm = norm(Omega - ADJ_x_in);
	 //////////////////////////////

	 logger.cleanup();
	 // Flip gradient
	 objgrad.reverseInPlace();
	 dDdy.reverseInPlace();
	 dDdtheta.reverseInPlace();

    dDdtheta(0) = dDdtheta(size1);
		dDdy(0) = dDdy(size1);

// logger.cleanup();
eta.reverseInPlace();     // remove?
zeta.reverseInPlace();   // remove?
dgdtheta_save.reverseInPlace();   // remove?

eta(size1) = eta(0);
zeta(size1) = zeta(0);
// dgdtheta_save(size1) = dgdtheta_save(0);

for (int n = 0; n < size3; ++n){

	dA(n) = 0.0;
	dB(n) = 0.0;
	dC(n) = 0.0;
	dD(n) = 0.0;
}

// switch between two types of Fourier modes
if (Fourier_mode > 0.5){
	// y = sum (An * sin + Bn * cos)
	// theta = sum (Cn * sin + Dn * cos)
	for (int n = 0; n < size3; ++n){
	  for (int k = 0; k < period; ++k){
	  /// add dgdA ...
			dA(n) = dA(n) + eta(k) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * dt;   /// keep dt or not
			dB(n) = dB(n) + eta(k) * cos(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * dt;
			dC(n) = dC(n) + zeta(k) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * dt;
			dD(n) = dD(n) + zeta(k) * cos(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * dt;

			dA_check(k) = dA(size3-1);
			increment(k) = eta(k) * sin(2 * PI * (2*size3) * (double)k*dt/((double)period*dt)) * dt;
			increment1(k) = eta(k) * sin(2 * PI * (2*size3) * (double)k*dt/((double)period*dt));
			// cout << " eta(k)... " << eta(k) << endl;
			// cout << " zeta(k)... " << zeta(k) << endl;
		}

	}

}
else if (Fourier_mode < 0.5){
	// y = sum (An * sin(t + Bn))
	// theta = sum (Cn * sin(t + Dn))

	for (int n = 0; n < size3; ++n){
	  for (int k = 0; k < period; ++k){
	  /// add dgdA ...
				dA(n) = dA(n) + eta(k) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt) + B(n)) * dt;   /// keep dt or not
				dB(n) = dB(n) + eta(k) * A(n) * cos(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt) + B(n)) * dt;
				dC(n) = dC(n) + zeta(k) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt) + D(n)) * dt;
				dD(n) = dD(n) + zeta(k) * C(n) * cos(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt) + D(n)) * dt;

				dA_check(k) = dA(size3-1);
				increment(k) = eta(k) * sin(2 * PI * (2*size3) * (double)k*dt/((double)period*dt) + B(size3-1)) * dt;
				increment1(k) = eta(k) * sin(2 * PI * (2*size3) * (double)k*dt/((double)period*dt) + B(size3-1));
				// cout << " eta(k)... " << eta(k) << endl;
				// cout << " zeta(k)... " << zeta(k) << endl;

		}

	}

}
// compute gradient with respect to T
/////
// eta.reverseInPlace();     // remove?
// zeta.reverseInPlace();   // remove?
// eta(0) = eta(size1);
// zeta(0) = zeta(size1);
/////
double dT1 = 0.0;
double dT2 = 0.0;
double sum1;
double sum2;
if (Fourier_mode > 0.5){
	// y = sum (An * sin + Bn * cos)
	// theta = sum (Cn * sin + Dn * cos)
	for (int k = 0; k <= period; ++k){
		sum1 = 0.0;
		for (int n = 0; n < size3; ++n){
			sum1 = sum1 + (A(n) * cos(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.)\
			 - B(n) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.));
			 // cout << " sum1 ... " << sum1 << endl;
			 // cout << " n ... " << n << endl;
		}
		// dT1 = dT1 + eta(k) * sum1 * dt;
		dT1 = dT1 + eta(k) * sum1 * dt;
		// cout << " k ... " << k << endl;
		// cout << " eta(k)... " << eta(k) << endl;
		// cout << " dT1 ... " << dT1 << endl;
	}

	for (int k = 0; k <= period; ++k){
		sum2 = 0.0;
		for (int n = 0; n < size3; ++n){
			sum2 = sum2 + (C(n) * cos(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.)\
			 - D(n) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.));
			 // cout << " sum2 ... " << sum2 << endl;
			 // cout << " n ... " << n << endl;
		}
		dT2 = dT2 + zeta(k) * sum2 * dt;
		// dT2 = dT2 + zeta(k) * sum2 * dt + dgdtheta_save(k) * sum2 * dt;
		// cout << " k ... " << k << endl;
		// cout << " (double)k*dt ... " << (double)k*dt << endl;
		// cout << " zeta(k)... " << zeta(k) << endl;
		// cout << " dT2 ... " << dT2 << endl;
	}

	// dT = (dT1 + dT2) + dgdT + g_T;  // add g(T) at line 1257 and 1269
	// dT = (dT1 + dT2) + dgdT + g_T;  // add g(T) at line 1257 and 1269
	 cout << " dT ... " << dT << " dT1 ... " << dT1 << " dT2 ... " << dT2 << " dgdT ... " << dgdT << " g_T ... " << g_T<< endl;
	// cout << " (double)period*dt ... " << (double)period*dt << endl;
	// cout << " pow((double)period*dt,-2.) ... " << pow((double)period*dt,-2.) << endl;
}

if (Fourier_mode < 0.5){
	// y = sum (An * sin(t + Bn))
	// theta = sum (Cn * sin(t + Dn))
	for (int k = 0; k <= period; ++k){
		sum1 = 0.0;
		for (int n = 0; n < size3; ++n){
			sum1 = sum1 + A(n) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt) + B(n)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.);
      // cout << " sum1 ... " << sum1 << endl;
			// cout << " n ... " << n << endl;
		}
		dT1 = dT1 + eta(k) * sum1 * dt;
		// cout << " k ... " << k << endl;
		// cout << " eta(k)... " << eta(k) << endl;
		// cout << " dT1 ... " << dT1 << endl;
	}

	for (int k = 0; k <= period; ++k){
		sum2 = 0.0;
		for (int n = 0; n < size3; ++n){
			sum2 = sum2 + C(n) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt) + D(n)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.);
      // cout << " sum2 ... " << sum2 << endl;
			// cout << " n ... " << n << endl;
		}
		dT2 = dT2 + zeta(k) * sum2 * dt;
		// cout << " k ... " << k << endl;
		// cout << " (double)k*dt ... " << (double)k*dt << endl;
		// cout << " zeta(k)... " << zeta(k) << endl;
		// cout << " dT2 ... " << dT2 << endl;
	}
  	dT = (dT1 + dT2) + dgdT + g_T;  // add g(T) at line 1257 and 1269
	// dT = (dT1 + dT2) + dgdT + g_T;
	cout << " dT ... " << dT << " dT1 ... " << dT1 << " dT2 ... " << dT2 << " dgdT ... " << dgdT << " g_T ... " << g_T<< endl;
	// cout << " (double)period*dt ... " << (double)period*dt << endl;
	// cout << " pow((double)period*dt,-2.) ... " << pow((double)period*dt,-2.) << endl;
}

/////
// eta.reverseInPlace();     // remove?
// zeta.reverseInPlace();   // remove?
/////
// create dDdy
// Save dDdt
ofstream myfile_dA_check("../examples/dA_check.txt");
if (myfile_dA_check.is_open())
{
// for (int count = 0; count < size; count++) {
for (int count = 0; count < size1; ++count) {
	myfile_dA_check << dA_check(count) << " ";
}
myfile_dA_check.close();
}
else cout << "Unable to open file";

// create dDdy
// Save dDdt
ofstream myfile_din("../examples/din_check.txt");
if (myfile_din.is_open())
{
// for (int count = 0; count < size; count++) {
for (int count = 0; count < size1; ++count) {
	myfile_din << increment(count) << " ";
}
myfile_din.close();
}
else cout << "Unable to open file";

// create dDdy
// Save dDdt
ofstream myfile_din1("../examples/din1_check.txt");
if (myfile_din1.is_open())
{
// for (int count = 0; count < size; count++) {
for (int count = 0; count < size1; ++count) {
	myfile_din1 << increment1(count) << " ";
}
myfile_din1.close();
}
else cout << "Unable to open file";


//
//// create dDdtheta
//std::ifstream file_dDdtheta("../examples/dDdtheta_pp.txt");
//double dddtheta[numSteps+1];
//// Update Error
//for (int i = 0; i <= numSteps; ++i)
// {
//	 file_dDdtheta >> dddtheta[i];
// }
//// Save dDdt
//ofstream myfile_dDdtheta("../examples/dDdtheta_pp.txt");
//if (myfile_dDdtheta.is_open())
//{
// // for (int count = 0; count < size; count++) {
// for (int count = 0; count <= numSteps; ++count) {
//	 myfile_dDdtheta << dDdtheta[count] << " ";
// }
// myfile_dDdtheta.close();
//}
//else cout << "Unable to open file";
//
//// create dDdy
//std::ifstream file_dDdy("../examples/dDdy_pp.txt");
//double dddy[numSteps+1];
//// Update Error
//for (int i = 0; i <= numSteps; ++i)
//	{
//		file_dDdy >> dddy[i];
//	}
//// Save dDdt
//ofstream myfile_dDdy("../examples/dDdy_pp.txt");
//if (myfile_dDdy.is_open())
//{
//	// for (int count = 0; count < size; count++) {
//	for (int count = 0; count <= numSteps; ++count) {
//		myfile_dDdy << dDdy[count] << " ";
//	}
//	myfile_dDdy.close();
//}
//else cout << "Unable to open file";




// Save dDdt
ofstream myfile_eta("../examples/" + outdir + "eta.txt");
if (myfile_eta.is_open())
{
 // for (int count = 0; count < size; count++) {
 for (int count = 0; count < numSteps; ++count) {
	 myfile_eta << eta(count) << " ";
 }
 myfile_eta.close();
}
else cout << "Unable to open file";

// Save dDdt
ofstream myfile_zeta("../examples/" + outdir + "zeta.txt");
if (myfile_zeta.is_open())
{
 // for (int count = 0; count < size; count++) {
 for (int count = 0; count < numSteps; ++count) {
	 myfile_zeta << zeta(count) << " ";
 }
 myfile_eta.close();
}
else cout << "Unable to open file";


delete solver;
return Omega - ADJ_x_in;

}

/////////////
// adjoint to compute the gradient

inline Vec Adj(const Mat& input) {

	int numCycles = Ini_step;
	double lift = 0.;
	double drag = 0.;
  Grid grid( nx, ny, ngrid, length, xOffset, yOffset);
	// Make a flat plate, length 1, with center at 1/4 chord
	RigidBody plate;
	plate.addLine(0, 0, 1, 0, grid.Dx());
	plate.setCenter(0, 0);

	//
	int TT = size1;
	int end_index = TT * numCycles;
 	// double pitchPeriod_temp = round(TT / dt);
 	// // TT = pitchPeriod_temp * dt;
	// int pitchPeriod = (int)pitchPeriod_temp;

	double pitch_in[size1*numCycles+1];
	double plunge_in[size1*numCycles+1];
	double motiontime_in[size1*numCycles+1];


	for(int i = 0; i<size1*numCycles+1; ++i){
				pitch_in[i] =  input(i,0);
				plunge_in[i] =  input(i,1);
				motiontime_in[i] = input(i,2);

				// pitch_in[i] =  MotionParameter(i,0);
				// plunge_in[i] =  MotionParameter(i,1);
				// motiontime_in[i] = MotionTime(i);


	}

	PitchPlungeT_mod motion(
		dt,
		motiontime_in,
		pitch_in,
		plunge_in,
		end_index
		// matrix.row(i)
		// matrix.col(j)

	);
	plate.setMotion(motion);
	Geometry geom;
	geom.addBody(plate);
	// geom.moveBodies(0);

	// Setup equations to solve

	double alpha = 0;  // angle of background flow
	int periodStart = 0;
	int period = TT;
	BaseFlow q_potential(grid, magnitude, alpha);

/////////////////////////////
	// Motion* mot = geom.transferMotion(); // pull motion from first RigidBody object
	// geom.transferCenter(xC, yC);  // pull center of motion from RigidBody object
	// q_potential.setMotion(*mot);
	// q_potential.setCenter(xC, yC);
//////////////////////////
	NavierStokesModel* model = NULL;
	AdjointIBSolver2* adjointsolver = NULL;
	IBSolver* solver = NULL;

  State x00(grid, geom.getNumPoints());
	//vector<State> x0(period, x00);
	x0_glob.resize(period + 1);
/////////////////
  Motion* mot = NULL; // this is for adjoint tests; get rid of later
///////////////
/////////////////////////////
	mot = geom.transferMotion(); // pull motion from first RigidBody object
	geom.transferCenter(xC, yC);  // pull center of motion from RigidBody object
	q_potential.setMotion(*mot);
	q_potential.setCenter(xC, yC);
//////////////////////////

//for (int i = 0; i < period; i++){
//	// sprintf("PBF", "PBF".c_str(), i + periodStart);
//	// if ( ! x0[i].load(pbffilename) ) {
//
//	cout << "i = " << std::to_string(i) << endl;
//
//	if (i < 10) {
//		cout << "first" << endl;
//		if ( ! x0[i].load(outdir_PBF + "pitch00" + std::to_string(i) + ".bin") ) {
//				cout << "base flow " << outdir_PBF + "pitch00" + std::to_string(i) + ".bin" << " failed to load. Exiting program." << endl;
//				exit(1);
//			 }
//
//	}
//	if(i >= 10 && i < 100){
//		cout << "second" << endl;
//		if ( ! x0[i].load(outdir_PBF + "pitch0" + std::to_string(i) + ".bin") ) {
//				cout << "base flow " << outdir_PBF + "pitch0" + std::to_string(i) + ".bin" << " failed to load. Exiting program." << endl;
//				exit(1);
//			 }
//
//	}
//	if(i >= 100){
//		cout << "third" << endl;
//		if ( ! x0[i].load(outdir_PBF + "pitch" + std::to_string(i) + ".bin") ) {
//				cout << "base flow " << outdir_PBF + "/pitch" + std::to_string(i) + ".bin" << " failed to load. Exiting program." << endl;
//				exit(1);
//			 }
//
//	}
//
//}


  //x00 = x0[0];
	for (int i = 0; i <= period; i++) {
		// sprintf("PBF", "PBF".c_str(), i + periodStart);
		// if ( ! x0[i].load(pbffilename) ) {

		cout << "i = " << std::to_string(i) << endl;

		if (i < 10) {
			cout << "first" << endl;
			if (!x0_glob[i].load(outdir_PBF + "pitch00" + std::to_string(i) + ".bin")) {
				cout << "base flow " << outdir_PBF + "pitch00" + std::to_string(i) + ".bin" << " failed to load. Exiting program." << endl;
				exit(1);
			}

		}
		if (i >= 10 && i < 100) {
			cout << "second" << endl;
			if (!x0_glob[i].load(outdir_PBF + "pitch0" + std::to_string(i) + ".bin")) {
				cout << "base flow " << outdir_PBF + "pitch0" + std::to_string(i) + ".bin" << " failed to load. Exiting program." << endl;
				exit(1);
			}

		}
		if (i >= 100) {
			cout << "third" << endl;
			if (!x0_glob[i].load(outdir_PBF + "pitch" + std::to_string(i) + ".bin")) {
				cout << "base flow " << outdir_PBF + "pitch" + std::to_string(i) + ".bin" << " failed to load. Exiting program." << endl;
				exit(1);
			}

		}
		// x0[i] = x00;
		// x0_glob.push_back(x0[i]);
   // x0_glob.push_back(x00);
		cout << "x0_glob size... " << x0_glob.size() << endl;
	}
	///////////////////////
	// mot = geom.transferMotion();
	////////////////////
	model = new NavierStokesModel( grid, geom, Reynolds);
	// model = new NavierStokesModel( grid, geom, Reynolds, q_potential );


	Scheme::SchemeType str2scheme( string integratorType );
  string integratorType = "rk3";
	// string integratorType = "AB2"; //////////////////////
	// string integratorType = "EULER"; //////////////////////
	Scheme::SchemeType schemeType = str2scheme( integratorType );

	adjointsolver = new AdjointIBSolver2( grid, *model, dt, schemeType, x0_glob, period, cost, CT, CP, *mot );
	// adjointsolver = new AdjointIBSolver2( grid, *model, dt, "rk3", x0, period, *mot );
	solver = adjointsolver;

	State x( grid, geom.getNumPoints() );
	x.omega = 0.;
	x.f = 0.;
	x.q = 0.;

	// if not the first ieration, load bin file that gives better initial gauess
	if (k_global != 0) {
		x.load(outdir_ADJ + "ADJ000.bin");
	}

	State Boundary(grid, geom.getNumPoints());
	Boundary.f = 0.;
	  // update the geometry to the current time
	// geom.moveBodies( x.time );

	// Initialize model and timestepper
	model->init();
	cout << "using " << solver->getName() << " timestepper" << endl;
	cout << "    dt = " << dt << "\n" << endl;

			solver->init();
			solver->save( "outdir_name" );

	// Calculate flux for state, in case only vorticity was saved
	if( ! q_potential.isStationary() ) {
			q_potential.setAlphaMag( x.time );
			alpha = q_potential.getAlpha();
			cout << "q_potential calculation" << endl;
	}
	// model->updateOperators( x.time ); /////////////////////////
	// model->refreshState( x );  ////////////////////////

	cout << endl << "Initial timestep = " << x.timestep << "\n" << endl;


	StateVector xv(x);

  int numSteps = numCycles * period;

	 // start adjoint solver
	 Vec objgrad;
	 int k;
	 objgrad.resize(numSteps+1); // vector to store gradient
	 double dx2 = grid.Dx() * grid.Dx();
	 double dx = grid.Dx();
	 double alpha_previous;
	 double beta_previous;
	 double alpha_adj;
	 double beta_adj;
	 double lift_adj = 0.;
	 double drag_adj = 0.;
	 double drag_adj_pre = 0.;
	 double y = 0.;
	 double y_pre = 0.;
	 double dgdT;
	 Vec dDdy(numSteps+1);
	 Vec dDdtheta(numSteps+1);
	 Vec eta_nonperiodic(size1 * num_cycles + 1);
	 Vec zeta_nonperiodic(size1 * num_cycles + 1);


   Vec dgdtheta_save(numSteps+1);


	 double dqdaN_s[numSteps+1];

	 double dqdbN_s[numSteps+1];

	 double dudaphi_s[numSteps+1];

	 double dudbphi_s[numSteps+1];


// tranform the flow from NS using unsteady base flow to steady base flow !!!!
//////////////////////
  for(int i = 0; i < numSteps; i++) {
	 // for(int i = 1; i <= numSteps; i++) {
		 if(i > period){
			 k = std::remainder(i,period);
		 }
		 if(k < 0){
			 k = k + period;
		 }


			 //k = period - xv.x.timestep - 1; // if not periodic
			 k = (period - (xv.x.timestep % period)) % period; // if periodic
			 double partial = 0.;
			 Flux dqdp(x0_glob[k].q);
			 TangentSE2 gg(0, 0, 0, 0, 1., 0);  // eq 14 in Daniel's note, add the one for eq 15 which is the second argument.
			 // (x, y, theta, xdot, ydot, thetadot).
			 dqdp.setFlow(gg, xC, yC);
			 double deriv1 = (-InnerProduct(dqdp, CrossProduct(xv.x.q, x0_glob[k].omega))) / dx2;
			 BoundaryVector dudp(xv.x.f.getNumPoints());
			 dudp = -model->flux2boundary(dqdp);
			 double deriv2 = InnerProduct(dudp, xv.x.f);

			 // logger.doOutput( q_potential, x);

			 objgrad(i) = partial + deriv1 + deriv2;


       ///////////////////
			 // a = hdot, b = thetadot

			 // y_pre = y;
			 // y = motion.gety(xv.x.time);
			 // double theta = motion.gettheta(xv.x.time);
			 // double thetadot = motion.getthetadot(xv.x.time);
			 // double ydot = motion.getydot(xv.x.time);
			 /////////////////
       y_pre = y;
			 y = motion.gety(k * dt);
			 double theta = motion.gettheta(k * dt);
			 double thetadot = motion.getthetadot(k * dt);
			 double ydot = motion.getydot(k * dt);
			 double y = motion.gety(k * dt);

/////////////////////////////
/////////////////////////////
////////////////////////////
///// dgd()
double dgdh = 0.;
double dgdhdot = 0.;
double dgdtheta = 0.;
double dgdthetadot = 0.;
dgdT = 0.;
double xF, yF, zM;
x0_glob[k].computeNetForce(xF, yF);
// to check if the time step is right or note can be removed
drag_adj_pre = drag_adj;
drag_adj = xF * cos(theta) - yF * sin(theta);
lift_adj = xF * sin(theta) + yF * cos(theta);

// compute the pitching moment zM
Boundary.f = plate.getPoints();
zM = 0.0;
for (int j = 0;j < geom.getNumPoints();++j) {
//zM = x.f(0, j) * grid.Dx() * grid.Dx() * (Boundary.f(1, j) - yC) * sin(theta) + x.f(1, j) * grid.Dx() * grid.Dx() * (Boundary.f(0, j) - xC) * cos(theta) + zM;
// zM = x.f(1, j) * grid.Dx() * grid.Dx() * (Boundary.f(0, j) - xC) + zM;
zM = x0_glob[k].f(1, j) * grid.Dx() * grid.Dx() * (Boundary.f(0, j) - xC) + zM;
// cout << " x0[k].f(1, j) ... " << x0[k].f(1, j) << endl;
// cout << " grid.Dx() ... " << grid.Dx() << endl;
// cout << " (Boundary.f(0, j) - xC) ... " << (Boundary.f(0, j) - xC) << endl;
// cout << " zM ... " << zM << endl;
}

if (cost == 0){
	// cost = 0, thrust
 // dgdtheta = -(-xF * sin(theta) - yF * cos(theta))*2; // * 2
 dgdtheta = (-xF * sin(theta) - yF * cos(theta)) * 2 / ((double)period*dt) * magnitude; // * 2  add /(period * dt) for everything related to g
 // dgdtheta = (xF * sin(-theta) + yF * cos(-theta))*2; // * 2
 dgdtheta_save(i) = dgdtheta;

 dgdthetadot = 0 / ((double)period*dt);
 dgdh = 0 / ((double)period*dt);
 dgdhdot = 0 / ((double)period*dt);
 // dgdT = (-1) * pow((double)period*dt,-2.) * CT * 2;
 dgdT = (-1) * CT / ((double)period * dt);
 // cout << "(double)period ... " << (double)period << endl;
 // cout << "(double)period*dt ... " << (double)period*dt << endl;
 // cout << "pow((double)period*dt,-2) ... " << pow((double)period*dt,-2) << endl;
 // cout << "CT ... " << CT << endl;
}

else if (cost == 1){
	// cost = 1, efficiency
	dgdtheta = 1/CP * (- 2 * sin(theta) *  xF - 2 * cos(theta) * yF) / ((double)period * dt) - CT/(pow(CP,2)) * (2 * sin(theta) * yF - 2 * cos(theta) * xF) * ydot / ((double)period * dt);
	// dgdtheta = (2 * sin(theta) * yF * ydot - 2 * cos(theta) * xF * ydot - 2 * sin(theta) * xF - 2 * cos(theta) * yF) / ((double)period*dt);
	// dgdtheta = (-2 * sin(theta) * yF * ydot + 2 * cos(theta) * xF * ydot - 2 * sin(theta) * xF - 2 * cos(theta) * yF);


	// dgdthetadot = -(-zM);  // -
	dgdthetadot = -(-(CT/(pow(CP,2)) * (-2 * zM) / ((double)period * dt)));
	// dgdthetadot = -(-2 * zM) / ((double)period*dt);  // -
	cout << " zM ... " << zM << endl;

	dgdh = 0 / ((double)period*dt);

	dgdhdot = -(-(CT/(pow(CP,2)) * (-2 *cos(theta) * yF - 2 * sin(theta) * xF) / ((double)period * dt)));
	// dgdhdot = -(-(CT/(pow(CP,2)) * (-2 *cos(theta) * yF - sin(theta) * xF)));
	// dgdhdot = -(-2 * cos(theta) * yF -2 * sin(theta) * xF) / ((double)period*dt);  /// -

	dgdT = 0;
	// dgdT = (-1) * pow((double)period*dt,-2) * (CP + CT) * dt * 2;


 // // Note that the dgd() in the note is -dgd()
 // dgdtheta = (-1) * ((xF * sin(theta) + yF * cos(theta)) * (-cos(theta) * yF * ydot -sin(theta) * xF * ydot - zM * thetadot)/(pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period\
 // 					- (-cos(theta) * xF + yF * sin(theta)) * (sin(theta) * yF * ydot -cos(theta) * xF * ydot)/(pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period);
 //
 // dgdthetadot = (-1) * (0/(pow(-cos(theta) * yF *ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period\
 // 						 - (-xF * cos(theta) + yF * sin(theta) * (-1) * zM)/(pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period);
 //
 // dgdh = 0;
 //
 // dgdhdot	= (-1) * (0/(pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period\
 // 					 - (-xF * cos(theta) + yF * sin(theta)) * (-cos(theta) * yF - sin(theta) * xF) / (pow(-cos(theta) * yF * ydot - sin(theta) * xF * ydot - zM * thetadot,2)) / period);
 //
 // dgdT = (-1) * pow(period,-2) * xF/(yF * ydot + zM * thetadot);
}
//////////////////
////////////////////////////
//////////////////////////
////////////////////////

//////////////////////////////
// alpha
State x_dqda(grid, geom.getNumPoints());
Flux dqda(x_dqda.q);

Flux::index ind;
for (int lev=0; lev<dqda.Ngrid(); ++lev) {

		double xdot_dqda = -sin(theta);
		double ydot_dqda = -cos(theta);
		// double xdot = 0;
		// double ydot = 1;
		double xdiff, ydiff;


		// double dx = dqda.Dx(lev);
		for(ind = dqda.begin(X); ind != dqda.end(X); ++ind) {
				xdiff = dqda.x(lev,ind) - xC;
				ydiff = dqda.y(lev,ind) - yC;
				// dqda(lev,ind) = (xdot -thetadot*ydiff)*dx;
				dqda(lev,ind) = (xdot_dqda)*dx;
		}
		for(ind = dqda.begin(Y); ind != dqda.end(Y); ++ind) {
				xdiff = dqda.x(lev,ind) - xC;
				ydiff = dqda.y(lev,ind) - yC;
				// dqda(lev,ind) = (ydot + thetadot*xdiff)*dx;
				dqda(lev,ind) = (ydot_dqda)*dx;
		}
}

// double deriv1_alpha = (-InnerProduct(dqda, CrossProduct(xv.x.q, x0[k].omega))) / dx2;  // cout this term and the one in next line
double deriv1_alpha = (InnerProduct(dqda, CrossProduct(xv.x.q, x0_glob[k].omega))) / dx2;  // cout this term and the one in next line ///////mod
BoundaryVector duda(xv.x.f.getNumPoints());
duda = -model->flux2boundary(dqda); // remove -
// xv.x.f = -1 * xv.x.f;   ////////////////////////////
double deriv2_alpha = -InnerProduct(duda, xv.x.f);
alpha_previous = alpha_adj;
alpha_adj = deriv1_alpha + deriv2_alpha;
double dalphadt = (alpha_adj - alpha_previous) / dt;  /// add -  !!!
////////////////////////////////
////////
// cout << " deriv1_alpha ... " << deriv1_alpha << endl;
// cout << "deriv2_alpha ... " << deriv2_alpha << endl;
/////////////////
State x_dqdb(grid, geom.getNumPoints());
Flux dqdb(x_dqdb.q);
Flux::index ind_1;
for (int lev=0; lev<dqdb.Ngrid(); ++lev) {

		double xdiff, ydiff;

		// double dx = dqda.Dx(lev);
		for(ind_1 = dqdb.begin(X); ind_1 != dqdb.end(X); ++ind_1) {
				xdiff = dqdb.x(lev,ind_1) - xC;
				ydiff = dqdb.y(lev,ind_1) - yC;
				// dqda(lev,ind) = (xdot -thetadot*ydiff)*dx;
				dqdb(lev,ind_1) = -(ydiff)*dx;
		}
		for(ind_1 = dqdb.begin(Y); ind_1 != dqdb.end(Y); ++ind_1) {
				xdiff = dqdb.x(lev,ind_1) - xC;
				ydiff = dqdb.y(lev,ind_1) - yC;
				// dqda(lev,ind) = (ydot + thetadot*xdiff)*dx;
				dqdb(lev,ind_1) = -(-xdiff)*dx;
		}
}
//////////////
       // Flux dqdb(x0[k].q);
       // TangentSE2 gg_b(0, 0, 0, 0, 0, -1.);  // eq 14 in Daniel's note, add the one for eq 15 which is the second argument. The right one
			 // // TangentSE2 gg_b(0, 0, 0, 0, 0, 1);  // eq 14 in Daniel's note, add the one for eq 15 which is the second argument. The right one
       // // (x, y, theta, xdot, ydot, thetadot).
			 // dqdb.setFlow(gg_b, xC, yC);
/////////////
			 // double deriv1_beta = (-InnerProduct(dqdb, CrossProduct(xv.x.q, x0[k].omega))) / dx2;
			 double deriv1_beta = (InnerProduct(dqdb, CrossProduct(xv.x.q, x0_glob[k].omega))) / dx2;  //////////////////////// mod
			 BoundaryVector dudb(xv.x.f.getNumPoints());
			 dudb = -model->flux2boundary(dqdb);
			 double deriv2_beta = -InnerProduct(dudb, xv.x.f);
       beta_previous = beta_adj;

			 beta_adj = deriv1_beta + deriv2_beta;
			 double dbetadt = (beta_adj - beta_previous) / dt;
			 /////////////////
			 // cout << " dalpha/dt ..." << dalphadt << " dbeta/dt ..." << dbetadt << endl;

///////////
State x_dqdh(grid, geom.getNumPoints());
Flux dqdh(x_dqdb.q);
Flux::index ind_2;
for (int lev=0; lev<dqdh.Ngrid(); ++lev) {

		double xdiff, ydiff;

		// double dx = dqda.Dx(lev);
		for(ind_2 = dqdh.begin(X); ind_2 != dqdh.end(X); ++ind_2) {
				xdiff = dqdh.x(lev,ind_2) - xC;
				ydiff = dqdh.y(lev,ind_2) - yC;

				dqdh(lev,ind_2) = -thetadot * sin(theta) * dx * 0;  //////////////////////
		}
		for(ind_2 = dqdh.begin(Y); ind_2 != dqdh.end(Y); ++ind_2) {
				xdiff = dqdh.x(lev,ind_2) - xC;
				ydiff = dqdh.y(lev,ind_2) - yC;

				dqdh(lev,ind_2) = thetadot * cos(theta) * dx * 0;  //////////////////
		}
}

BoundaryVector dudh(xv.x.f.getNumPoints());
dudh = -model->flux2boundary(dqdh);
///////////
			 // Flux dqdh(x0[k].q);
			 // dqdh = dqdh * 0;
			 // BoundaryVector dudh(xv.x.f.getNumPoints());
			 // dudh = -model->flux2boundary(dqdh);
////////////
/////////////////////
/////////////////
State x_dqdtheta(grid, geom.getNumPoints());
Flux dqdtheta(x_dqdtheta.q);
Flux::index ind_3;
for (int lev=0; lev<dqdtheta.Ngrid(); ++lev) {

		double xdiff, ydiff;

		// double dx = dqda.Dx(lev);
		for(ind_3 = dqdtheta.begin(X); ind_3 != dqdtheta.end(X); ++ind_3) {
				xdiff = dqdtheta.x(lev,ind_3) - xC;
				ydiff = dqdtheta.y(lev,ind_3) - yC;
				// dqdtheta(lev,ind_3) = (-sin(theta) - thetadot*xdiff  -ydot*cos(theta))*dx;
				dqdtheta(lev,ind_3) = (-sin(theta)-ydot*cos(theta))*dx;  ////////////////////////
		}
		for(ind_3 = dqdtheta.begin(Y); ind_3 != dqdtheta.end(Y); ++ind_3) {
				xdiff = dqdtheta.x(lev,ind_3) - xC;
				ydiff = dqdtheta.y(lev,ind_3) - yC;
				// dqdtheta(lev,ind_3) = (-cos(theta) + thetadot*ydiff +ydot*sin(theta))*dx;
				dqdtheta(lev,ind_3) = (-cos(theta)+ydot*sin(theta))*dx;   /////////////////////////
		}
}
/////////////////////
       // Flux dqdtheta(x0[k].q);
			 // TangentSE2 gg_theta(0, 0, 0, -sin(theta)-ydot*cos(theta) , -cos(theta)+ydot*sin(theta), 0);  // eq 14 in Daniel's note, add the one for eq 15 which is the second argument. The right one
       // // (x, y, theta, xdot, ydot, thetadot).
			 // dqdtheta.setFlow(gg_theta, xC, yC);
///////////////////
			 BoundaryVector dudtheta(xv.x.f.getNumPoints());
			 dudtheta = -model->flux2boundary(dqdtheta);

			 double delta_h;
			 double delta_theta;
			 // delta_h = dgdh-(InnerProduct(dqdh, CrossProduct(xv.x.q, x0[k].omega)))/dx2-InnerProduct(dudh, xv.x.f)-dalphadt; //
			 // x0[k].omega = x0[k].omega + 2 * thetadot;///////////////////////////////////////////////////////
			 delta_h = dgdh+(InnerProduct(dqdh, CrossProduct(xv.x.q, x0_glob[k].omega)))/dx2-InnerProduct(dudh, xv.x.f)-dalphadt; // ////////////// mod
			 // delta_theta = dgdtheta-(InnerProduct(dqdtheta, CrossProduct( xv.x.q, x0[k].omega)))/dx2-InnerProduct(dudtheta, xv.x.f)-dbetadt;
			 delta_theta = dgdtheta+(InnerProduct(dqdtheta, CrossProduct( xv.x.q, x0_glob[k].omega)))/dx2-InnerProduct(dudtheta, xv.x.f)-dbetadt; ////////////// mod
			 cout << " delta_h ..." << delta_h << " delta_theta ..." << delta_theta << endl;
			 //////////
			 // cout << " dgdh ... " << dgdh << endl;
			 // double dqdaN = (InnerProduct(dqdh, CrossProduct(xv.x.q, x0[k].omega)))/dx2;
			 // cout << " dqda*N ... " << dqdaN << endl;
			 //
			 // double dqdbN = (InnerProduct(dqdtheta, CrossProduct(xv.x.q, x0[k].omega)))/dx2;
			 // cout << " dqdb*N ... " << dqdbN << endl;
			 //
			 // double dudaphi = InnerProduct(dudh, xv.x.f);
			 // cout << " duda*phi ... " << dudaphi << endl;
			 // double dudbphi = InnerProduct(dudtheta, xv.x.f);
			 // cout << " dudb*phi ... " << dudbphi << endl;
			 // cout << " dalphadt ... " << dalphadt << endl;
			 // double Ssin = sin(xv.x.time/2*(2*PI));
			 // cout << " sin(t) ... " << Ssin << endl;

       //////////
       dDdy(i) = delta_h;
       dDdtheta(i) = delta_theta;
			 // dqdaN_s[i] = dqdaN;
			 // dqdbN_s[i] = dqdbN;
			 // dudaphi_s[i] = dudaphi;
			 // dudbphi_s[i] = dudbphi;
			 if (i <= period) {
				 eta_nonperiodic(i) = delta_h;
				 zeta_nonperiodic(i) = delta_theta;
			 }
			 cout << " dDdy ..." << dDdy(i) << endl;
			 solver->advance(xv.x);

	 }

	 // Flip gradient
	 objgrad.reverseInPlace();
	 dDdy.reverseInPlace();
	 dDdtheta.reverseInPlace();
	 eta_nonperiodic.reverseInPlace();
	 zeta_nonperiodic.reverseInPlace();
  ///////////////////////
	///////////////////////
	//////////////////////
	/////
	double dT1 = 0.0;
	double dT2 = 0.0;
	double sum1;
	double sum2;
	if (Fourier_mode > 0.5){
		// y = sum (An * sin + Bn * cos)
		// theta = sum (Cn * sin + Dn * cos)
		for (int k = 0; k <= period; ++k){
			sum1 = 0.0;
			for (int n = 0; n < size3; ++n){
				sum1 = sum1 + (A(n) * cos(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.)\
				 - B(n) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.));

			}
			// dT1 = dT1 + eta(k) * sum1 * dt;
			dT1 = dT1 + eta_nonperiodic(k) * sum1 * dt;

		}

		for (int k = 0; k <= period; ++k){
			sum2 = 0.0;
			for (int n = 0; n < size3; ++n){
				sum2 = sum2 + (C(n) * cos(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.)\
				 - D(n) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.));

			}
			dT2 = dT2 + zeta_nonperiodic(k) * sum2 * dt;
			// dT2 = dT2 + zeta(k) * sum2 * dt + dgdtheta_save(k) * sum2 * dt;

		}

		dT_nonperiodic = (dT1 + dT2) + dgdT + g_T;  // add g(T) at line 1257 and 1269
		// dT = (dT1 + dT2) + dgdT + g_T;  // add g(T) at line 1257 and 1269
		// cout << " dT ... " << dT << " dT1 ... " << dT1 << " dT2 ... " << dT2 << " dgdT ... " << dgdT << " g_T ... " << g_T<< endl;
		cout << " dT_nonperiodic ... " << dT_nonperiodic << " dT1 ... " << dT1 << " dT2 ... " << dT2 << " dgdT ... " << dgdT << " g_T ... " << g_T << endl;
		// cout << " (double)period*dt ... " << (double)period*dt << endl;
		// cout << " pow((double)period*dt,-2.) ... " << pow((double)period*dt,-2.) << endl;
	}

	if (Fourier_mode < 0.5){
		// y = sum (An * sin(t + Bn))
		// theta = sum (Cn * sin(t + Dn))
		for (int k = 0; k <= period; ++k){
			sum1 = 0.0;
			for (int n = 0; n < size3; ++n){
				sum1 = sum1 + A(n) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt) + B(n)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.);
			}
			dT1 = dT1 + eta_nonperiodic(k) * sum1 * dt;

		}

		for (int k = 0; k <= period; ++k){
			sum2 = 0.0;
			for (int n = 0; n < size3; ++n){
				sum2 = sum2 + C(n) * sin(2 * PI * (2*n+1) * (double)k*dt/((double)period*dt) + D(n)) * 2 * PI * (2*n+1) * (double)k*dt * (-1) * pow((double)period*dt,-2.);

			}
			dT2 = dT2 + zeta_nonperiodic(k) * sum2 * dt;

		}
	  	dT_nonperiodic = (dT1 + dT2) + dgdT + g_T;  // add g(T) at line 1257 and 1269
		// dT = (dT1 + dT2) + dgdT + g_T;
		cout << " dT_nonperiodic ... " << dT_nonperiodic << " dT1 ... " << dT1 << " dT2 ... " << dT2 << " dgdT ... " << dgdT << " g_T ... " << g_T<< endl;
		// cout << " (double)period*dt ... " << (double)period*dt << endl;
		// cout << " pow((double)period*dt,-2.) ... " << pow((double)period*dt,-2.) << endl;
	}

	//////////////////////
	///////////////////////
	/////////////////////

	 ///////////////////
	 // Convert Vec variables to State
 	int M = 0;
 	//Vec Omega_in(((nx-1) * (ny-1) * 3) * ngrid + geom.getNumPoints() * 2);
	Vec Omega_in(((nx - 1)* (ny - 1) + (nx + 1) * ny + nx * (ny + 1))* ngrid + geom.getNumPoints() * 2);

 	for (int k = 0; k < ngrid; ++k){
 		for (int i = 1; i < nx; ++i) {
 			for (int j = 1; j < ny; ++j) {
 				Omega_in(M) = xv.x.omega(k, i, j);
 				M = M + 1;
 			}

 		}

 		for (int i = 0; i <= nx; ++i) {
 			for (int j = 0; j < ny; ++j) {
 				Omega_in(M) = xv.x.q(k, 0, i, j);
 				M = M + 1;
 			}

 		}

 		for (int i = 0; i < nx; ++i) {
 			for (int j = 0; j <= ny; ++j) {
 				Omega_in(M) = xv.x.q(k, 1, i, j);
 				M = M + 1;
 			}

 		}

 	}


 	for (int i = 0; i < geom.getNumPoints(); ++i) {
 		Omega_in(M) = xv.x.f(0, i);
 		M = M + 1;
 	}

 	for (int i = 0; i < geom.getNumPoints(); ++i) {
 		Omega_in(M) = xv.x.f(1, i);
 		M = M + 1;
 	}

	linear_solver_ns::Gmres<Vec> gmres(innerproduct, 5, 1.e-3);   //1.e-3 1.e-10
	// constexpr double tol{ Newton_tol };
	constexpr double tol{ 1.e-8 };   //1.e-3 1.e-6  1.e-8
	constexpr int max_iter{ 15 };
	constexpr double jacobian_dx{ 1.e-3 };  //1.e-3 1.e-6
	std::unique_ptr<newton_ns::Newton<Vec>> newton;
		newton.reset(
			new newton_ns::Newton<Vec>{ ADJ_F, gmres, norm, tol, max_iter,
					jacobian_dx, 1, true }
		);

	//newton.reset(
	//	new newton_ns::NewtonArmijo<Vec>{ F, gmres, norm, tol, max_iter,
	//			jacobian_dx, 1, 10, 0.1, 0.5, 1.e-4, true }
	//);

	try {
		newton->solve(Omega_in);
	}

	catch (newton_ns::NewtonError& ne) {
		cerr << "Newton error: " << ne.what() << endl;
		exit(2);
	}
	catch (std::exception& e) {
		cerr << "Other exception: " << e.what() << endl;
		exit(3);
	}
	 ///////////////////

	 //////////////
	   save_flag = 1;
	   ADJ_F(Omega_in);
	 	save_flag = 0;
	 /////////////

// logger.cleanup();

delete solver;
return objgrad;
// return gradient;

}
//

class Drag_average
 {
 private:
	  int n; // the number of 
	  int id; // id of processor
	  int N; // number of processors
 public:
	 Drag_average(int n_, int id_, int N_) : n(n_), id(id_), N(N_) {}
	 double operator()(const Vec& x, Vec& grad) // x is the control parameters, grad is the gradient
	 {


		Vec tau(2);
		double grad_norm = 10; // use this to save drag or efficiency for forward run 
		double CT_plus; // use this to save the perturbed effeciency or drag for forward run
		double CT_0; // use this to save the nonperturbed effecincy or drag for forward run
		double CT_temp; // use this to combine and save both perturbed and nonperturbed drag or effeciency
    for(int i = 0; i < size3; ++i){
			A(i) = x(i);
			B(i) = x(i + size3);
			C(i) = x(i + size3 * 2);
			D(i) = x(i + size3 * 3);
			dt = x(n-1);
		}

	/////////////////////////
		 ///////////////
	if (id == 0) {

		// create A
		std::ifstream file_Aa("../examples/" + outdir + "A" + std::to_string(k_global) + ".txt");
		double Aa[size3];
		// Update Error
		for (int i = 0; i < size3; ++i)
		{
			file_Aa >> Aa[i];
		}
		// Save A
		ofstream myfile_Aa("../examples/" + outdir + "A" + std::to_string(k_global) + ".txt");
		if (myfile_Aa.is_open())
		{
			// for (int count = 0; count < size; count++) {
			for (int count = 0; count < size3; ++count) {
				myfile_Aa << A(count) << " ";
			}
			myfile_Aa.close();
		}
		else cout << "Unable to open file";

		// create B
		std::ifstream file_Ba("../examples/" + outdir + "B" + std::to_string(k_global) + ".txt");
		double Ba[size3];
		// Update Error
		for (int i = 0; i < size3; ++i)
		{
			file_Ba >> Ba[i];
		}
		// Save B
		ofstream myfile_Ba("../examples/" + outdir + "B" + std::to_string(k_global) + ".txt");
		if (myfile_Ba.is_open())
		{
			// for (int count = 0; count < size; count++) {
			for (int count = 0; count < size3; ++count) {
				myfile_Ba << B(count) << " ";
			}
			myfile_Ba.close();
		}
		else cout << "Unable to open file";


		// create C
		std::ifstream file_Ca("../examples/" + outdir + "C" + std::to_string(k_global) + ".txt");
		double Ca[size3];
		// Update Error
		for (int i = 0; i < size3; ++i)
		{
			file_Ca >> Ca[i];
		}
		// Save C
		ofstream myfile_Ca("../examples/" + outdir + "C" + std::to_string(k_global) + ".txt");
		if (myfile_Ca.is_open())
		{
			// for (int count = 0; count < size; count++) {
			for (int count = 0; count < size3; ++count) {
				myfile_Ca << C(count) << " ";
			}
			myfile_Ca.close();
		}
		else cout << "Unable to open file";


		// create D
		std::ifstream file_Da("../examples/" + outdir + "D" + std::to_string(k_global) + ".txt");
		double Da[size3];
		// Update Error
		for (int i = 0; i < size3; ++i)
		{
			file_Da >> Da[i];
		}
		// Save D
		ofstream myfile_Da("../examples/" + outdir + "D" + std::to_string(k_global) + ".txt");
		if (myfile_Da.is_open())
		{
			// for (int count = 0; count < size; count++) {
			for (int count = 0; count < size3; ++count) {
				myfile_Da << D(count) << " ";
			}
			myfile_Da.close();
		}
		else cout << "Unable to open file";

		// create T
		std::ifstream file_Ta("../examples/" + outdir + "T" + std::to_string(k_global) + ".txt");
		double Ta;
		// Update Error

		file_Ta >> Ta;

		// Save D
		ofstream myfile_Ta("../examples/" + outdir + "T" + std::to_string(k_global) + ".txt");
		if (myfile_Ta.is_open()) {
			myfile_Ta << dt << " ";
			myfile_Ta.close();
		}
		else cout << "Unable to open file";
	}
	/////////////////
		   // compute CT, CP with tau = tau + dt for FD of dJdtau
		 	if (id == 1){
		 		T_perturb = 1;
		 		dt = dt + delta_T;
		 		tau(0) = (double)size1 * dt;
		 		tau(1) = (double)size1 * dt;

		 		MotionParameterTPerturb.resize(size1 * Ini_step + 1,size2);
		 		MotionTimeTPerturb.resize(size1 * Ini_step + 1);
		 		MotionParameterTimeTPerturb.resize(size1 * Ini_step + 1,size2+1);
		 	// generate time series motion data with Fourier modes
		 	MotionParameterTimeTPerturb = MotionGenerator( A, B, C, D, Fourier_mode, tau, size1*Ini_step, size2, dt, PI);
		 	MotionParameterTPerturb.col(0) = MotionParameterTimeTPerturb.col(0);
		 	MotionParameterTPerturb.col(1) = MotionParameterTimeTPerturb.col(1);
		 	MotionTimeTPerturb = MotionParameterTimeTPerturb.col(2);

		 	outdir_PBF = outdir_PBF_global + "T_perturb/";
		 		CT_plus = drag_int_fun(MotionParameterTimeTPerturb);

		 	  // CT_plus = 1;
		 		cout << "CT_plus ... " << CT_plus << endl;
		 		cout << "CT ..." << CT << endl;
		 	}

		 if (id == 0){

		 	// compute base flow
		 	T_perturb = 0;
		 	dt = dt;
		 	tau(0) = (double)size1 * dt;
		 	tau(1) = (double)size1 * dt;

		 	MotionParameter.resize(size1 * Ini_step + 1,size2);
		 	MotionTime.resize(size1 * Ini_step + 1);
		 	MotionParameterTime.resize(size1 * Ini_step + 1,size2+1);

		 // generate time series motion data with Fourier modes
		 MotionParameterTime = MotionGenerator( A, B, C, D, Fourier_mode, tau, size1*Ini_step, size2, dt, PI);
		 MotionParameter.col(0) = MotionParameterTime.col(0);
		 MotionParameter.col(1) = MotionParameterTime.col(1);
		 MotionTime = MotionParameterTime.col(2);
		 /////////////////
				// create dDdy
		 std::ifstream file_ya("../examples/" + outdir + "y" + std::to_string(k_global) + ".txt");
		 double ya[size1 * Ini_step + 1];
		 // Update Error
		 for (int i = 0; i <= size1 * Ini_step; ++i)
		 {
			 file_ya >> ya[i];
		 }
		 // Save dDdt
		 ofstream myfile_ya("../examples/" + outdir + "y" + std::to_string(k_global) + ".txt");
		 if (myfile_ya.is_open())
		 {
			 // for (int count = 0; count < size; count++) {
			 for (int count = 0; count < MotionParameter.rows(); ++count) {
				 myfile_ya << MotionParameter(count, 1) << " ";
			 }
			 myfile_ya.close();
		 }
		 else cout << "Unable to open file";


		 // create dDdy
		 std::ifstream file_thetaa("../examples/" + outdir + "theta" + std::to_string(k_global) + ".txt");
		 double thetaa[size1 * Ini_step + 1];
		 // Update Error
		 for (int i = 0; i <= size1 * Ini_step; ++i)
		 {
			 file_thetaa >> thetaa[i];
		 }
		 // Save dDdt
		 ofstream myfile_thetaa("../examples/" + outdir + "theta" + std::to_string(k_global) + ".txt");
		 if (myfile_thetaa.is_open())
		 {
			 // for (int count = 0; count < size; count++) {
			 for (int count = 0; count < MotionParameter.rows(); ++count) {
				 myfile_thetaa << MotionParameter(count, 0) << " ";
			 }
			 myfile_thetaa.close();
		 }
		 else cout << "Unable to open file";

		 //////////////////


		 outdir_PBF = outdir_PBF_global;
		 	CT_0 = drag_int_fun(MotionParameterTime);
		 	// CT_0 = 0;
		 	cout << "CT_0 ... " << CT_0 << endl;
		 	cout << "CT ..." << CT << endl;


			// create dDdy
			std::ifstream file_La("../examples/" + outdir + "L" + std::to_string(k_global) + ".txt");
			double La[size1];
			// Update Error
			for (int i = 0; i < size1; ++i)
			{
				file_La >> La[i];
			}
			// Save dDdt
			ofstream myfile_La("../examples/" + outdir + "L" + std::to_string(k_global) + ".txt");
			if (myfile_La.is_open())
			{
				// for (int count = 0; count < size; count++) {
				// for (int count = 0; count <= size1 * 8; ++count) {
				for (int count = 0; count < size1; ++count) {
					myfile_La << L_save(count) << " ";
				}
				myfile_La.close();
			}
			else cout << "Unable to open file";
			///
			// create dDdy
			std::ifstream file_Draga("../examples/" + outdir + "Drag" + std::to_string(k_global) + ".txt");
			double Draga[size1];
			// Update Error
			for (int i = 0; i < size1; ++i)
			{
				file_Draga >> Draga[i];
			}
			// Save dDdt
			ofstream myfile_Draga("../examples/" + outdir + "Drag" + std::to_string(k_global) + ".txt");
			if (myfile_Draga.is_open())
			{
				// for (int count = 0; count < size; count++) {
				// for (int count = 0; count <= size1 * 8; ++count) {
				for (int count = 0; count < size1; ++count) {
					myfile_Draga << D_save(count) << " ";
				}
				myfile_Draga.close();
			}
			else cout << "Unable to open file";
			///
				 // create dDdy
			std::ifstream file_za("../examples/" + outdir + "zm" + std::to_string(k_global) + ".txt");
			double za[size1];
			// Update Error
			for (int i = 0; i < size1; ++i)
			{
				file_za >> za[i];
			}
			// Save dDdt
			ofstream myfile_za("../examples/" + outdir + "zm" + std::to_string(k_global) + ".txt");
			if (myfile_za.is_open())
			{
				// for (int count = 0; count < size; count++) {
				// for (int count = 0; count <= size1 * 8; ++count) {
				for (int count = 0; count < size1; ++count) {
					myfile_za << z_m_save(count) << " ";
				}
				myfile_za.close();
			}
			else cout << "Unable to open file";


		 // 	// Use adjoint to compute the gradient for the rest parameters
		 	Adj(MotionParameterTime);
		 // // compute dJdT with finite differnece

		 }


		 if (id == 0) {
		 	CT_temp = CT_0; // assign drag or effeiciency from forward run to processor 0
		 }
		 if (id == 1) {
		 	CT_temp = CT_plus; // assign drag or effeciency from forward run to processor 1
		 }
		 cout << "Line 2979" << endl;

		 /////////////

		 cout << "Line 2990" << endl;
		 double* rbuf_temp;
		 rbuf_temp = (double*)malloc(N * 1 * sizeof(double));

		 MPI_Barrier(MPI_COMM_WORLD);
		 MPI_Allgather(&CT_temp, 1, MPI_DOUBLE, rbuf_temp, 1, MPI_DOUBLE, MPI_COMM_WORLD); // gether all CT_temp to processor 0 and assign to rbuf_temp[0] and rbuf_temp[1] to processor 0
		 cout << "Line 2995" << endl;

		  //dT = (rbuf_temp[1] - rbuf_temp[0]) / delta_T;
		  dT = (rbuf_temp[1] - rbuf_temp[0]) / delta_T * dt;
		  ////////////////////
		  // in case dt has too many digits, truncate it and remove some digits
		  dT = trunc(dT * 1.e7)/(1.e7);
		  ////////////////////
		  cout << "f(t+dt)-f ..." << (rbuf_temp[1] - rbuf_temp[0]) << endl;
		  cout << "dt ..." << dt;
		 //dT = (rbuf_temp[1] - rbuf_temp[0]) / delta_T/(size1*dt)/(size1);
		 // dT = dT / dt;
		 // dT = dT/(abs(dT)) * ceil(abs(dT/dt));
		 //

		 ///////////////////////////

		  if (id == 0){
		 	if (cost == 0){
		 		grad_norm = CT; // drag from forward run
		 	}

		 	else if (cost == 1){
		 		grad_norm = CT/CP; // effeciency from forward run
		 	}

		  }
//////
		  MPI_Barrier(MPI_COMM_WORLD);
		  // broadcast everyting from processor 0 to all processors
		 MPI_Bcast(A.data(), A.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		 MPI_Bcast(B.data(), B.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		 MPI_Bcast(C.data(), C.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		 MPI_Bcast(D.data(), D.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		 MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		 MPI_Bcast(dA.data(), dA.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		 MPI_Bcast(dB.data(), dB.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		 MPI_Bcast(dC.data(), dC.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		 MPI_Bcast(dD.data(), dD.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		 MPI_Bcast(&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		 MPI_Bcast(&grad_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		 /////////////

		meandrag[k_global] = grad_norm; // save drag or effeciency for each gradient iteration
    cout << " k_global ... " << k_global << endl;
		////////////////

		 cout << "Line 3187" << endl;
		 /////////////////
		 ////////////////

		 MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////////////

        if (id == 0)
		{
 			// create drag
 			std::ifstream file_Drag("../examples/" + outdir + "Drag_opt.txt");
 			double D_temp[num_iter * num_iter_linsearch];
 			// Update Error
 			for (int i = 0; i < num_iter * num_iter_linsearch; ++i)
 				{
 					file_Drag >> D_temp[i];
 				}
 			// Save drag
 			ofstream myfile_Drag("../examples/" + outdir + "Drag_opt.txt");
 			if (myfile_Drag.is_open())
 			{
 				// for (int count = 0; count < size; count++) {
 				for (int count = 0; count < num_iter * num_iter_linsearch; ++count) {
 					myfile_Drag << meandrag[count] << " ";
 				}
 				myfile_Drag.close();
 			}
 			else cout << "Unable to open file";
 			///////////////////

 		 // create dA
 		 std::ifstream file_dAa("../examples/" + outdir + "dA" + std::to_string(k_global) + ".txt");
 		 double dAa[size3];
 		 // Update Error
 		 for (int i = 0; i < size3; ++i)
 			 {
 				 file_dAa >> dAa[i];
 			 }
 		 // Save dA
 		 ofstream myfile_dAa("../examples/" + outdir + "dA" + std::to_string(k_global) + ".txt");
 		 if (myfile_dAa.is_open())
 		 {
 			 // for (int count = 0; count < size; count++) {
 			 for (int count = 0; count < size3; ++count) {
 				 myfile_dAa << dA(count) << " ";
 			 }
 			 myfile_dAa.close();
 		 }
 		 else cout << "Unable to open file";

 		 // create dB
 		 std::ifstream file_dBa("../examples/" + outdir + "dB" + std::to_string(k_global) + ".txt");
 		 double dBa[size3];
 		 // Update Error
 		 for (int i = 0; i < size3; ++i)
 			 {
 				 file_dBa >> dBa[i];
 			 }
 		 // Save dB
 		 ofstream myfile_dBa("../examples/" + outdir + "dB" + std::to_string(k_global) + ".txt");
 		 if (myfile_dBa.is_open())
 		 {
 			 // for (int count = 0; count < size; count++) {
 			 for (int count = 0; count < size3; ++count) {
 				 myfile_dBa << dB(count) << " ";
 			 }
 			 myfile_dBa.close();
 		 }
 		 else cout << "Unable to open file";


 		 // create dC
 		 std::ifstream file_dCa("../examples/" + outdir + "dC" + std::to_string(k_global) + ".txt");
 		 double dCa[size3];
 		 // Update Error
 		 for (int i = 0; i < size3; ++i)
 			 {
 				 file_dCa >> dCa[i];
 			 }
 		 // Save dC
 		 ofstream myfile_dCa("../examples/" + outdir + "dC" + std::to_string(k_global) + ".txt");
 		 if (myfile_dCa.is_open())
 		 {
 			 // for (int count = 0; count < size; count++) {
 			 for (int count = 0; count < size3; ++count) {
 				 myfile_dCa << dC(count) << " ";
 			 }
 			 myfile_dCa.close();
 		 }
 		 else cout << "Unable to open file";


 		 // create dD
 		 std::ifstream file_dDa("../examples/" + outdir + "dD" + std::to_string(k_global) + ".txt");
 		 double dDa[size3];
 		 // Update Error
 		 for (int i = 0; i < size3; ++i)
 			 {
 				 file_dDa >> dDa[i];
 			 }
 		 // Save dD
 		 ofstream myfile_dDa("../examples/" + outdir + "dD" + std::to_string(k_global) + ".txt");
 		 if (myfile_dDa.is_open())
 		 {
 			 // for (int count = 0; count < size; count++) {
 			 for (int count = 0; count < size3; ++count) {
 				 myfile_dDa << dD(count) << " ";
 			 }
 			 myfile_dDa.close();
 		 }
 		 else cout << "Unable to open file";


 		 // create dT
 		 std::ifstream file_dTa("../examples/" + outdir + "dT" + std::to_string(k_global) + ".txt");
 		 double dTa;
 		 // Update Error

 			 file_dTa >> dTa;

 		 // Save dD
 		 ofstream myfile_dTa("../examples/" + outdir + "dT" + std::to_string(k_global) + ".txt");
 		 if (myfile_dTa.is_open())
 		 {
 			 myfile_dTa << dT << " ";
 			 myfile_dTa.close();
 		 }
 		 else cout << "Unable to open file";
 		 //////////////////////////////
 		}

		if (id == 1){


			// create drag
			std::ifstream file_Drag("../examples/" + outdir + "Drag_opt_id1.txt");
			double D_temp[num_iter * num_iter_linsearch];
			// Update Error
			for (int i = 0; i < num_iter * num_iter_linsearch; ++i)
				{
					file_Drag >> D_temp[i];
				}
			// Save drag
			ofstream myfile_Drag("../examples/" + outdir + "Drag_opt_id1.txt");
			if (myfile_Drag.is_open())
			{
				// for (int count = 0; count < size; count++) {
				for (int count = 0; count < num_iter * num_iter_linsearch; ++count) {
					myfile_Drag << meandrag[count] << " ";
				}
				myfile_Drag.close();
			}
			else cout << "Unable to open file";
			///////////////////

		 //////////////////////////////
		}



		 // }
		 // cout << "A0 ... " << A(0) << endl;
		 // cout << "A1 ... " << A(1) << endl;
		 // cout << "A2 ... " << A(2) << endl;
		 // cout << "A3 ... " << A(3) << endl;
		 // cout << "A4 ... " << A(4) << endl;

     k_global = k_global + 1;
		 double fx = 0.0;
     for (int i = 0; i < size3; ++i){
			 grad(i) = dA(i);

			 grad(i+size3) = dB(i);

			 grad(i+size3*2) = dC(i);

			 grad(i+size3*3) = dD(i);
		 }
		 grad(n-1) = dT;
     fx = grad_norm;
     MPI_Barrier(MPI_COMM_WORLD);
		 return fx;
	 }
 };
//////////////////////
int main(int argc, char* argv[]) {
// initialize the parameters
initialization();
  Vec tau(2);

  Mat theta_temp(size1 * Ini_step + 1, size3);
  Mat y_temp(size1 * Ini_step + 1, size3);
	double grad_norm_pre;
  double grad_norm = 10;
	double CT_plus;
	double CT_0;
	double CT_temp;

	//

	double meandrag[num_iter * num_iter_linsearch];
	tau(0) = (double)size1 * dt; // period for pitching
	tau(1) = (double)size1 * dt; // period for heaving


	for (int i = 0; i < size3; ++i){
		A(i) = A_0(i);
		B(i) = B_0(i);
		C(i) = C_0(i);
		D(i) = D_0(i);
	}

// generate time series motion data with Fourier modes
MotionParameterTime = MotionGenerator( A, B, C, D, Fourier_mode, tau, size1*Ini_step, size2, dt, PI); // double check tau(0) and tau(1) in MotionGenerator, there seems a mistake!!
MotionParameter.col(0) = MotionParameterTime.col(0);
MotionParameter.col(1) = MotionParameterTime.col(1);
MotionTime = MotionParameterTime.col(2);

// make a folder to save all the parameters
AddSlashToPath(outdir);
mkdir(outdir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
AddSlashToPath(outdir_PBF_global);
AddSlashToPath(outdir_ADJ);
// Save command
ofstream myfile_command("../examples/" + outdir + "command.txt");
if (myfile_command.is_open())
{
	// for (int count = 0; count < size; count++) {

	myfile_command << "numner of total sample points, size1=" + std::to_string(size1) << "\n";
  myfile_command << "2 for pitch and plunge, size2=" + std::to_string(size2) << "\n";
	myfile_command << "number of Fourier modes, size3=" + std::to_string(size3) << "\n";
	myfile_command << "dt=" + std::to_string(dt) << "\n";
	myfile_command << "Reynolds number =" + std::to_string(Reynolds) << "\n";
	myfile_command << "Freestream velocity, magnitude =" + std::to_string(magnitude) << "\n";
	myfile_command << "number of points in horizontal direction, nx =" + std::to_string(nx) << "\n";
	myfile_command << "number of points in vertical direction, ny =" + std::to_string(ny) << "\n";
	myfile_command << "number of grids, ngrid =" + std::to_string(ngrid) << "\n";
	myfile_command << "save all grids or not, TecplotAllGrids =" + std::to_string(TecplotAllGrids) << "\n";
	myfile_command << "length of the smallest window, length =" + std::to_string(length) << "\n";
	myfile_command << "offset in horizontal direction, xOffset =" + std::to_string(xOffset) << "\n";
	myfile_command << "offset in vertical direction, yOffset =" + std::to_string(yOffset) << "\n";
	myfile_command << "length of each cell in vertical direction, dy =" + std::to_string(dy) << "\n";
	myfile_command << "pitching axis in horizontal direction, xC =" + std::to_string(xC) << "\n";
	myfile_command << "pitching axis in vertical direction, yC =" + std::to_string(yC) << "\n";
	myfile_command << "total upper limit for pitching amplitude, theta_up =" + std::to_string(theta_up) << "\n";
	myfile_command << "upper limit of the heaving Fourier mode A, A_lim =" + std::to_string(A_lim) << "\n";
	myfile_command << "upper limit of the heaving Fourier mode B, B_lim =" + std::to_string(B_lim) << "\n";
	myfile_command << "upper limit of the pitching Fourier mode C, C_lim =" + std::to_string(C_lim) << "\n";
	myfile_command << "upper limit of the pitching Fourier mode D, D_lim =" + std::to_string(D_lim) << "\n";
	myfile_command << "upper limit of the pitching Fourier mode D, dt_lim_low =" + std::to_string(dt_lim_low) << "\n";
	myfile_command << "learning rate, alpha_step =" + std::to_string(alpha_step) << "\n";
	myfile_command << "maximum number of iterations, num_iter =" + std::to_string(num_iter) << "\n";
	myfile_command << "initial condition of the heaving Fourier mode A, A_0 =" + std::to_string(A_0(0)) << "\n";
	myfile_command << "initial condition of the heaving Fourier mode B, B_0 =" + std::to_string(B_0(0)) << "\n";
	myfile_command << "initial condition of the pitching Fourier mode C, C_0 =" + std::to_string(C_0(0)) << "\n";
	myfile_command << "initial condition of the pitching Fourier mode D, D_0 =" + std::to_string(D_0(0)) << "\n";
	myfile_command << "Type of Fourier mode, 0 for Asin(t+b), 1 for Asin(t)+Bsin(t), Fourier_mode =" + std::to_string(Fourier_mode) << "\n";
	myfile_command << "Type of cost function, 0 for drag, 1 efficiency, cost =" + std::to_string(cost) << "\n";
	myfile_command.close();
}
else cout << "Unable to open file";

// initialize MPI
int id;
int ierr;
int p;
int N = 2;
ierr = MPI_Init(NULL, NULL);

//  Get the number of processes.
//
ierr = MPI_Comm_size(MPI_COMM_WORLD, &p);
//
//  Get the individual process ID.
//
ierr = MPI_Comm_rank(MPI_COMM_WORLD, &id);

const int n = size3 * 4 + 1;
LBFGSBParam<double> param;
param.epsilon = 1e-3;
param.max_iterations = num_iter;
param.max_linesearch = num_iter_linsearch;
//param.ftol = 1e-4;
//param.wolfe = 0.9;
LBFGSpp::LBFGSBSolver<double> solver(param);
Drag_average fun(n,id,N);

// Variable bounds
Vec lb(n);
Vec ub(n);

for (int i = 0; i < size3; ++i) {
	ub(i) = A_lim;
	ub(i + size3) = B_lim;
	ub(i + size3 * 2) = C_lim;
	ub(i + size3 * 3) = D_lim;
	ub(n - 1) = dt_lim_up;

	lb(i) = -A_lim;
	lb(i + size3) = -B_lim;
	lb(i + size3 * 2) = -C_lim;
	lb(i + size3 * 3) = -D_lim;
	lb(n - 1) = dt_lim_low;
}

Vec y(n);

// initialize
for (int i = 0; i < size3; ++i) {
	y(i) = A(i);

	y(i + size3) = B(i);

	y(i + size3 * 2) = C(i);

	y(i + size3 * 3) = D(i);
}
y(n - 1) = dt;

double fx;
MPI_Barrier(MPI_COMM_WORLD);
int niter = solver.minimize(fun, y, fx, lb, ub, outdir);
std::cout << niter << " iterations" << std::endl;
std::cout << "x = \n" << y.transpose() << std::endl;
std::cout << "f(x) = " << fx << std::endl;
MPI_Barrier(MPI_COMM_WORLD);
////////////////////
MPI_Finalize();
	return 0;
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
