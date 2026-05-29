#ifndef _IBSOLVER_H_
#define _IBSOLVER_H_

// header file for storing integration schemes

#include <string>
#include <vector>
#include "Scheme.h"
#include "Scalar.h"
#include "State.h"
#include "Grid.h"
#include "NavierStokesModel.h"

namespace ibpm{

class ProjectionSolver;
class Motion;


// Base class
class IBSolver{
public:
    IBSolver( const Grid& grid,
		NavierStokesModel& model,
		double dt,
		Scheme::SchemeType scheme
	);
    IBSolver( const Grid& grid,
             NavierStokesModel& model,
             double dt,
             Scheme::SchemeType scheme,
             double tol
             );
	virtual ~IBSolver();
    void init();
    void reset();
	bool load(const string& basename);
	bool save(const string& basename);
	string getName();
    double getTimestep();
	void advance( State& x );
	void advance( State& x, const Scalar& Bu );
	virtual void advanceSubstep( State& x, const Scalar& nonlinear, int i ); // virtual for SFD and adjoint2
    void setTol( double tol );

protected:
	// methods
	virtual Scalar N( const State& x ) const = 0;
	ProjectionSolver* createSolver(double beta);
	void createAllSolvers();
	void deleteAllSolvers();

	// data
	const Grid& _grid;
	Scheme _scheme;
	string _name;
	double _dt;
	NavierStokesModel& _model;
	Scalar _Nprev;
	Scalar _Ntemp;
	bool _oldSaved;
    vector < ProjectionSolver* > _solver;
    double _tol;
};

// =============== //
// Derived classes //
// =============== //

class NonlinearIBSolver : public IBSolver {
public:
	NonlinearIBSolver(
		Grid& grid,
		NavierStokesModel& model,
		double dt,
		Scheme::SchemeType scheme
        ) :
        IBSolver( grid, model, dt, scheme ) { };

    NonlinearIBSolver(
        Grid& grid,
        NavierStokesModel& model,
        double dt,
        Scheme::SchemeType scheme,
        double tol
        ) :
        IBSolver( grid, model, dt, scheme, tol ) { };

protected:
	Scalar N( const State& x ) const;
};

class LinearizedIBSolver : public IBSolver {
public:
	LinearizedIBSolver(
		Grid& grid,
		NavierStokesModel& model,
		double dt,
		Scheme::SchemeType scheme,
		const State& baseFlow) :
		IBSolver( grid, model, dt, scheme ),
		_x0( baseFlow ) { }

    LinearizedIBSolver(
        Grid& grid,
        NavierStokesModel& model,
        double dt,
        Scheme::SchemeType scheme,
        double tol,
        const State& baseFlow) :
        IBSolver( grid, model, dt, scheme, tol ),
        _x0( baseFlow ) { }

protected:
	Scalar N(const State& x) const;

private:
	State _x0;
};

class AdjointIBSolver : public IBSolver {
public:
	AdjointIBSolver(
		Grid& grid,
		NavierStokesModel& model,
		double dt,
		Scheme::SchemeType scheme,
		const State& baseFlow)  :
		IBSolver( grid, model, dt, scheme ),
		_x0( baseFlow ) { }

	AdjointIBSolver(
        Grid& grid,
        NavierStokesModel& model,
        double dt,
        Scheme::SchemeType scheme,
        double tol,
        const State& baseFlow)  :
        IBSolver( grid, model, dt, scheme, tol ),
        _x0( baseFlow ) { }

protected:
	Scalar N(const State& x) const;

private:
	State _x0;
};

class AdjointIBSolver2 : public IBSolver {
  public:
   AdjointIBSolver2(
      Grid& grid,
      NavierStokesModel& model,
      double dt,
      Scheme::SchemeType scheme,
      const vector<State> x0periodic,
      const int period,
			const int cost,
      const double CT,
      const double CP,
      const Motion& motion) :
   IBSolver( grid, model, dt, scheme),
      _x0periodic( x0periodic ),
      _period( period ),
			_cost( cost ),
      _CT( CT ),
      _CP( CP ){
         assert(_period == static_cast<int>(x0periodic.size()));
         _motion = motion.clone();
      }

   AdjointIBSolver2(
      Grid& grid,
      NavierStokesModel& model,
      double dt,
      Scheme::SchemeType scheme,
      double tol,
      const vector<State> x0periodic,
      const int period,
			const int cost,
      const double CT,
      const double CP,
      const Motion& motion) :
   IBSolver( grid, model, dt, scheme, tol ),
      _x0periodic( x0periodic ),
      _period( period ),
			_cost( cost ),
      _CT( CT ),
      _CP( CP ){
         assert(_period == static_cast<int>(x0periodic.size()));
         _motion = motion.clone();
      }

  protected:
   Scalar N(const State& x) const;
   void advanceSubstep( State& x, const Scalar& nonlinear, int i );

  private:
   const vector<State> _x0periodic;
   const int _period;
	 const int _cost;
   const double _CT;
   const double _CP;
   Motion *_motion;
};

//! Navier-Stokes equations linearized about a periodic orbit.
class LinearizedPeriodicIBSolver : public IBSolver {
public:
	LinearizedPeriodicIBSolver(
		Grid& grid,
		NavierStokesModel& model,
		double dt,
		Scheme::SchemeType scheme,
		const vector<State> x0periodic,
		const int period )  :
		IBSolver( grid, model, dt, scheme ),
		_x0periodic( x0periodic ),
		_period( period )
	{
		assert(_period == static_cast<int>(x0periodic.size()));
	}

	LinearizedPeriodicIBSolver(
        Grid& grid,
        NavierStokesModel& model,
        double dt,
        Scheme::SchemeType scheme,
        double tol,
        const vector<State> x0periodic,
        const int period )  :
    IBSolver( grid, model, dt, scheme, tol ),
    _x0periodic( x0periodic ),
    _period( period )
	{
		assert(_period == static_cast<int>(x0periodic.size()));
	}

protected:
	Scalar N(const State& x) const;

private:
	const vector<State> _x0periodic;
	const int _period;
};

//! Navier-Stokes equations linearized about a periodic orbit.
class SFDSolver : public IBSolver {
public:
	SFDSolver(
		Grid& grid,
		NavierStokesModel& model,
		double dt,
		Scheme::SchemeType scheme,
		double Delta,
		double chi ) :
		IBSolver( grid, model, dt, scheme ),
		_Delta( Delta ),
		_chi( chi ),
		_xhat( _grid, _model.getNumPoints() ),
		_omegaTemp( grid ),
		_rhsPrev( grid ),
		_xhatSaved( false ),
		_rhsSaved( false ) { }

    void saveFilteredState( string outdir, string name, string numDigitInFileName );
    void loadFilteredState( string icFile );

protected:
	Scalar N(const State& x) const;
	void advanceSubstep( State& x, const Scalar& nonlinear, int i );

private:
	double _Delta;			// inverse of cutoff frequency
	double _chi;			// sfd gain
	State _xhat;
	Scalar _omegaTemp;
	Scalar _rhsPrev;
	bool _xhatSaved;
	bool _rhsSaved;
};


} // ibpm

#endif
