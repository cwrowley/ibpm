// update state x by Newton with know periodic, run IBPM for a few cycles for IC
// Using bracketing golden slection method to minimize the drag
// The parameters of motion are packed into one array for future multi-dimension optimization
// change the frequency to period for int number of steps

// Princeton University
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$
#ifndef POSolver
#define POSolver
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
//#include <vector>
//#include <iterator>

#include "StateVector.h"
#include "Gmres.h"
#include "Newton.h"
#include "NewtonArmijo.h"
#include "exceptions.h"

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
//using namespace Eigen;
//using namespace newton_ns; //

using Vec = Eigen::VectorXd;
using Scal = ibpm::Scalar;
// Define innerproduct
inline double innerproduct(const Vec& x, const Vec& y) {
	return y.dot(x);
}
// Define norm
inline double norm(const Vec& x) {
	return x.norm();
}
// Function to compute pitching for one cycle
inline Vec F(const Vec& Omega_in) {
	cout << "Pitching flat plate example\n";

	// lift and drag integral of drag
	double lift = 0.;
	double drag = 0.;
	double drag_int = 0;

	// Setup grid
	int nx = 200;
	int ny = 200;
	int ngrid = 1;
	double length = 4.0;
	double xOffset = -1;
	double yOffset = -2;
	Grid grid(nx, ny, ngrid, length, xOffset, yOffset);
  // Time step
	double dt = 0.001;

	// Load parameter for the pithcing plunging motion
	// para[pitching amplitude, pitching period, pitching phase, plunging amplitude, plunging frequency, plunging phase]
	std::ifstream file_para("../examples/para.txt");
	double para[6];
	int size_p = 6;

	for (int i = 0; i < size_p; ++i)
	{
		file_para >> para[i];
	}
    double A = para[0];
	double T = para[1];
	cout << "T..." << T << endl;

	// Convert T to freq

	double freq = 1 / (T * dt);

	// Make a flat plate, length 1, with center at 1/4 chord
	RigidBody plate;
	plate.addLine(0, 0, 1, 0, grid.Dx());
	plate.setCenter(0, 0);

	// Set the motion
	double amplitude = A;
	//double freq = 1;
	double phase = 0;
    //bool ubf = true;
	PitchPlunge motion(amplitude, freq, phase, 0, 0, 0);
	plate.setMotion(motion);
	Geometry geom;
	geom.addBody(plate);
	//geom.moveBodies(0);

	// Setup equations to solve
	double Reynolds = 100;
	double magnitude = 1;
	double alpha = 0;  // angle of background flow
	BaseFlow q_potential(grid, magnitude, alpha);
	double xC = 0;
	double yC = 0;

	Motion* m = geom.transferMotion(); // pull motion from first RigidBody object
	geom.transferCenter(xC, yC);  // pull center of motion from RigidBody object
	q_potential.setMotion(*m);
	q_potential.setCenter(xC, yC);

	cout << "Setting up Navier Stokes model..." << flush;
	NavierStokesModel model(grid, geom, Reynolds, q_potential);
	model.init();
	cout << "done" << endl;
	cout << "amplitude = " << A << endl;
	cout << "amplitude = " << amplitude << endl;

	// Setup timestepper
	NonlinearIBSolver solver(grid, model, dt, Scheme::AB2);
	solver.init();

	// Build the state variable, zero initial conditions
	State x(grid, geom.getNumPoints());
	x.omega = 0.;
	x.q = 0.;
	x.f = 0.;
	// Convert Vec variables to State
	int N = 0;
	for (int k = 0; k < ngrid;++k) {

		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				x.omega(k, i, j) = Omega_in(N);   //x.omega(0,i,j)
				N = N + 1;
			}

		}

		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				x.q(k, 0, i, j) = Omega_in(N);   //x.omega(0,i,j)
				N = N + 1;
			}

		}

		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
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

	cout << "Omega Newton" << Omega_in(2) << endl;
	// Define Fx to save the final value of State x
	State Fx(grid, geom.getNumPoints());
	Fx.omega = 0.;
	Fx.q = 0.;
	Fx.f = 0.;

	// Create output directory, if does not already exist
	mkdir("PO_test_out_11_ic_time_memo_UBF_IB_freq", S_IRWXU | S_IRWXG | S_IRWXO);

	// Setup output routines
	//OutputForce force("PO_test_out_11_ic_time_memo_UBF_IB_freq/force.dat");
	//OutputTecplot tecplot("PO_test_out_11_ic_time_memo_UBF_IB_freq/pitch%03d.plt", "Pitching plate, step %03d", 1);
	OutputRestart restart("PO_test_out_11_ic_time_memo_UBF_IB_freq/pitch%03d.bin");
	Logger logger;
	//Output Tecplot file every few timesteps
	//logger.addOutput(&tecplot, 1);
	//logger.addOutput(&force, 1);
	logger.addOutput(&restart, 100);
	logger.init();
	logger.doOutput(x);


	// Step
	const double PI = 4. * atan(1.);

/*	double cycle = 1 / freq / dt; */
	int numSteps = T;    //+1?
	for (int i = 1; i <= numSteps; ++i) {
		double theta = amplitude * sin(2 * PI * freq * x.time);
		cout << "step " << setw(4) << i
			<< "  time = " << setw(5) << x.time
			<< "  theta = " << theta << endl;

		solver.advance(x);
		Fx.omega = x.omega;
		Fx.q = x.q;
		Fx.f = x.f;
		double xF, yF;
		x.computeNetForce(xF, yF);

		q_potential.setAlphaMag(x.time);
		alpha = q_potential.getAlpha();
		drag = xF * cos(alpha) + yF * sin(alpha);
		lift = xF * -1. * sin(alpha) + yF * cos(alpha);

		// Compute drag integral
		drag_int = drag_int + drag;

		cout << " x force : " << setw(16) << drag * 2 << " , y force : "
			<< setw(16) << lift * 2 << "\n";

		logger.doOutput( q_potential, x);
	}

	// Convert State variables in to Vec
	int M = 0;
	Vec Omega(((nx-1) * (ny-1) * 3) * ngrid + geom.getNumPoints() * 2);
	for (int k = 0;k < ngrid; ++k) {
		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				Omega(M) = Fx.omega(k, i, j);
				M = M + 1;
			}

		}

		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				Omega(M) = Fx.q(k, 0, i, j);
				M = M + 1;
			}

		}

		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				Omega(M) = Fx.q(k, 1, i, j);
				M = M + 1;
			}

		}
	}
	

	for (int i = 0; i < geom.getNumPoints(); ++i) {
		Omega(M) = Fx.f(0, i);
		M = M + 1;
	}

	for (int i = 0; i < geom.getNumPoints(); ++i) {
		Omega(M) = Fx.f(1, i);
		M = M + 1;
	}

	//for (int i = 0; i < _numPoints; ++i) {
	//	x.f(0, i) = Omega_in(N);   //x.omega(0,i,j)
	//	N = N + 1;
	//}

	drag_int = drag_int / T;
	logger.cleanup();
	// Save drag integral
	ofstream file_drag_int;
	file_drag_int.open("drag_int.txt");
	file_drag_int << drag_int;
	file_drag_int.close();
	cout << "drag_in: " << drag_int << endl;


	// Compute the norm of x_final-x_initial
	double f_norm = norm(Omega - Omega_in);

	// Load step
	std::ifstream input_file_step("../examples/step.txt");
	int step;
	input_file_step >> step;
	// Load Error
	std::ifstream file_Err("../examples/ErrOut.txt");
	double Err[800];
	int size = 800;
	// Update Error
	for (int i = 0; i < size; ++i)
		{
			file_Err >> Err[i];
		}

	// Save Error
	Err[step] = { f_norm };
	ofstream myfile("../examples/ErrOut.txt");
	if (myfile.is_open())
	{
		// for (int count = 0; count < size; count++) {
		for (int count = 0; count < size; ++count) {
			myfile << Err[count] << " ";
		}
		myfile.close();
	}
	else cout << "Unable to open file";
	cout << "Error" << Err[step] << endl;

	step = step + 1;
	// Save step
	ofstream file;
	file.open("step.txt");
	file << step;
	file.close();
	cout << "Step_in: " << step << endl;


	return Omega - Omega_in;
	//return Omega;

}
// Newton and drag integral
inline double darg_int_fun(double T) {

	using namespace std::chrono;
	// Start timing
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	// Start Newton solver after a few IBPM cycles
	int Ini_step = 8;
	//
	using std::cout;
	using std::cerr;
	using std::endl;
	using std::string;
	using std::ofstream;

	// lift and drag
	double lift = 0.;
	double drag = 0.;

	// Time step
	double dt = 0.001;
	// Convert T to freq
	double freq = 1 / (T * dt);
	// Setup grid
	int nx = 200;
	int ny = 200;
	int ngrid = 1;
	double length = 4.0;
	double xOffset = -1;
	double yOffset = -2;
	Grid grid(nx, ny, ngrid, length, xOffset, yOffset);
	// Make a flat plate, length 1, with center at 1/4 chord
	RigidBody plate;
	plate.addLine(0, 0, 1, 0, grid.Dx());
	plate.setCenter(0, 0);

	// Set the motion
	double amplitude = 0.25;
	//double freq = 1;
	double phase = 0;
	PitchPlunge motion(amplitude, freq, phase, 0, 0, 0);
	// Save motion configuration
	const int size_p = 6;
	double para[6] = {amplitude, T, phase, 0, 0, 0};

	ofstream myfile_para;
	myfile_para.open("../examples/para.txt");

	//for (int count = 0; count < size_p; count++) {
	for (int count = 0; count < size_p; ++count) {
		myfile_para << para[count] << " ";
	}
	myfile_para.close();

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
	double drag_int = 0;
	ofstream file_drag_int;
	file_drag_int.open("../examples/drag_int.txt");
	file_drag_int << drag_int;
	file_drag_int.close();


	// Set plate motion
	plate.setMotion(motion);
	Geometry geom;
	geom.addBody(plate);
	//geom.moveBodies(0);

	// Setup equations to solve
	double Reynolds = 100;
	double magnitude = 1;
	double alpha = 0;  // angle of background flow
	BaseFlow q_potential(grid, magnitude, alpha);
	double xC = 0;
	double yC = 0;
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
	cout << "amplitude = " << amplitude << endl;

	// Setup timestepper
	NonlinearIBSolver solver(grid, model, dt, Scheme::AB2);
	solver.init();
	//
	State x(grid, geom.getNumPoints());
	x.omega = 0.1;
	x.q = 0.1;
	x.f = 0.1;
	//

	const double PI = 4. * atan(1.);

	//double cycle = 1 / freq / dt;
	//int cycle_int = (int)cycle;
	int numSteps = T * Ini_step;
	for (int i = 1; i <= numSteps; ++i) {
		double theta = amplitude * sin(2 * PI * freq * x.time);
		cout << "step " << setw(4) << i
			<< "  time = " << setw(5) << x.time
			<< "  theta = " << theta << endl;

		solver.advance(x);
		double xF, yF;
		x.computeNetForce(xF, yF);

		q_potential.setAlphaMag(x.time);
		alpha = q_potential.getAlpha();
		drag = xF * cos(alpha) + yF * sin(alpha);
		lift = xF * -1. * sin(alpha) + yF * cos(alpha);

	}

	//
	//string icfile = "pitch200.bin";
	//x.load(icfile);
	//cout << "Loading initial condition from file: " << icfile << endl;
	//

	// Convert Vec variables to State
	int M = 0;
	Vec Omega_in(((nx-1) * (ny-1) * 3) * ngrid + geom.getNumPoints() * 2);
	for (int k = 0;k < ngrid; ++k) {
		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				Omega_in(M) = x.omega(k, i, j);
				M = M + 1;
			}

		}

		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				Omega_in(M) = x.q(k, 0, i, j);
				M = M + 1;
			}

		}

		for (int i = 1; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
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

	// Newton-GMRES solver
	linear_solver_ns::Gmres<Vec> gmres(innerproduct, 5, 1.e-3);
	constexpr double tol{ 1.e-4 };
	constexpr int max_iter{ 50 };
	constexpr double jacobian_dx{ 1.e-3 };
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
	std::ifstream input_file_drag_int("../examples/drag_int.txt");
	//double drag_int;
	input_file_drag_int >> drag_int;

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	// End timing
	std::cout << "It took me " << time_span.count() << " seconds.";
	std::cout << std::endl;

	cout << "drag_int: " << drag_int << endl;
	return drag_int;
}
#endif // POSolver
