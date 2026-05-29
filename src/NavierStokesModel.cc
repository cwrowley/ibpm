// NavierStokesModel.cc
//
// Description: Implementation of the NavierStokesModel class
//
// Author(s):
// Clancy Rowley
//
// Date: 3 Jul 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "NavierStokesModel.h"

namespace ibpm {

	NavierStokesModel::NavierStokesModel(
	const Grid& grid,
	const Geometry& geometry,
	double Reynolds,
	const BaseFlow& q_potential
	) :
    _grid( grid ),
    _geometry( geometry ),
    _regularizer( grid, geometry ),
    _baseFlow( q_potential ),
    _ReynoldsNumber( Reynolds ),
    _poisson( grid ),
    _hasBeenInitialized( false )
	{}

    NavierStokesModel::NavierStokesModel(
    const Grid& grid,
    const Geometry& geometry,
    double Reynolds
    ) :
    _grid( grid ),
    _geometry( geometry ),
    _regularizer( grid, geometry ),
    _baseFlow( grid ),
    _ReynoldsNumber( Reynolds ),
    _poisson( grid ),
    _hasBeenInitialized( false )
    {
        _baseFlow.setFlux(0.);
    }

    NavierStokesModel::~NavierStokesModel() {}

    void NavierStokesModel::init() {
        if ( _hasBeenInitialized ) return;  // do only once
        // Update regularizer
        _regularizer.update();
        _hasBeenInitialized = true;
    }

    bool NavierStokesModel::isTimeDependent() const {
        bool flag = false;
        if( (!_geometry.isStationary()) || (!_baseFlow.isStationary()) ) flag = true;
        return flag;
    }

    bool NavierStokesModel::geTimeDependent() const {
        bool flag = false;
        if( (!_geometry.isStationary()) ) flag = true;
        return flag;
    }

    bool NavierStokesModel::bfTimeDependent() const {
        bool flag = false;
        if( (!_baseFlow.isStationary()) ) flag = true;
        return flag;
    }

    int NavierStokesModel::getNumPoints() const {
        return _geometry.getNumPoints();
    }

    double NavierStokesModel::getAlpha() const {
        return 1. / _ReynoldsNumber;
    }

    // Return the boundary velocities minus the base flow velocity at the boundary
    BoundaryVector NavierStokesModel::getConstraints() const {
        BoundaryVector b = _geometry.getVelocities();
        BoundaryVector b0 = getBaseFlowBoundaryVelocities();
        b -= b0;
        return b;
    }

    BoundaryVector NavierStokesModel::getPoints() const {
			BoundaryVector coords = _geometry.getPoints();
			return coords;
		}

    void NavierStokesModel::updateOperators( double time ) {
        if( bfTimeDependent() ) _baseFlow.moveFlow(time);
        if( geTimeDependent() ) {
            _geometry.moveBodies(time);
            _regularizer.update();
        }
    }

    void NavierStokesModel::B(const BoundaryVector& f, Scalar& omega ) const {
        assert( _hasBeenInitialized );
        Flux q = _regularizer.toFlux( f );
        Curl( q, omega );
    }

    void NavierStokesModel::C(const Scalar& omega, BoundaryVector& f) const {
        assert( _hasBeenInitialized );
		Flux q(_grid);
		computeFluxWithoutBaseFlow( omega, q );
		f = _regularizer.toBoundary( q );
	}

	void NavierStokesModel::computeFluxWithoutBaseFlow(const Scalar& omega,
													   Flux& q ) const {
		assert( _hasBeenInitialized );
		Scalar streamfunction = vorticityToStreamfunction( omega );
		Curl( streamfunction, q );
	}

	void NavierStokesModel::computeFlux(const Scalar& omega, Flux& q ) const {
		assert( _hasBeenInitialized );
		computeFluxWithoutBaseFlow( omega, q );
		q += _baseFlow.getFlux();
	}

	void NavierStokesModel::refreshState( State& x ) const {
		computeFlux( x.omega, x.q );
	}

	// Convert vorticity omega into streamfunction psi:
	//    Laplacian psi = - omega
	Scalar NavierStokesModel::vorticityToStreamfunction( const Scalar& omega ) const {
		assert( _hasBeenInitialized );
		// Solve L psi = omega, with zero Dirichlet bc's
		Scalar psi = -1. * omega;
		psi.coarsify();
		_poisson.solve( psi, psi );
		return psi;
	}

	BoundaryVector NavierStokesModel::getBaseFlowBoundaryVelocities() const {
		assert( _hasBeenInitialized );
		BoundaryVector velocity = _regularizer.toBoundary( _baseFlow.getFlux() );
		return velocity;
	}

	Flux NavierStokesModel::force2flux(BoundaryVector &f) {
		 Flux q = _regularizer.toFlux( f );
		 return q;
	}

	BoundaryVector NavierStokesModel::flux2boundary(Flux &q) {
		 BoundaryVector f = _regularizer.toBoundary( q );
		 return f;
	}


} // namespace ibpm
