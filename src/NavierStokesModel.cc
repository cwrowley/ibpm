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
    const Flux& q_potential
    ) :
    _grid( grid ),
    _geometry( geometry ),
    _regularizer( grid, geometry ),
    _baseFlow( q_potential ),
    _ReynoldsNumber( Reynolds ),
    _poisson( grid.Nx(), grid.Ny(), grid.Dx() ),
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
    _poisson( grid.Nx(), grid.Ny(), grid.Dx() ),
    _hasBeenInitialized( false )
    {
    _baseFlow = 0.;
}

NavierStokesModel::~NavierStokesModel() {}
    
void NavierStokesModel::init() {
    if ( _hasBeenInitialized ) return;  // do only once
    // Update regularizer
    _regularizer.update();    
    _hasBeenInitialized = true;
}

bool NavierStokesModel::isTimeDependent() const {
    return ! _geometry.isStationary();
}

int NavierStokesModel::getNumPoints() const {
    return _geometry.getNumPoints();
}
 
// Return the boundary velocities minus the base flow velocity at the boundary
BoundaryVector NavierStokesModel::getConstraints() const {
    BoundaryVector b = _geometry.getVelocities();
    BoundaryVector b0 = getBaseFlowBoundaryVelocities();
    b -= b0;
    return b;
}
    
void NavierStokesModel::updateOperators( double time ) {
    _geometry.moveBodies(time);
    _regularizer.update();
}
    
double NavierStokesModel::getAlpha() const {
    return 1. / _ReynoldsNumber;
}
    
void NavierStokesModel::B(const BoundaryVector& f, Scalar& gamma ) const {
    assert( _hasBeenInitialized );
    Flux q = _regularizer.toFlux( f );
    gamma = Curl( q );
}

void NavierStokesModel::C(const Scalar& gamma, BoundaryVector& f) const {
    assert( _hasBeenInitialized );
    Flux q(_grid);
    computeFluxWithoutBaseFlow( gamma, q );
    f = _regularizer.toBoundary( q );
}

void NavierStokesModel::computeFluxWithoutBaseFlow(const Scalar& gamma,
                                                   Flux& q ) const {
    assert( _hasBeenInitialized );
    Scalar streamfunction = gammaToStreamfunction( gamma );
    q = Curl( streamfunction );
}

void NavierStokesModel::computeFlux(const Scalar& gamma, Flux& q ) const {
    assert( _hasBeenInitialized );
    computeFluxWithoutBaseFlow( gamma, q );
    q += _baseFlow;
}

void NavierStokesModel::refreshState( State& x ) const {
    computeFlux( x.gamma, x.q );
}

// Convert circulation gamma into streamfunction psi:
//    Laplacian psi = - omega
//    gamma = omega * dx^2
//    psi = -1 / dx^2 * Laplacian^{-1} gamma    
Scalar NavierStokesModel::gammaToStreamfunction(const Scalar& gamma) const {
    assert( _hasBeenInitialized );
    Scalar psi( gamma.getGrid() );
    // Solve L psi = gamma, with zero Dirichlet bc's
    _poisson.solve( gamma, psi );
    psi *= -1. / ( gamma.Dx() * gamma.Dx() );
    return psi;
}

BoundaryVector NavierStokesModel::getBaseFlowBoundaryVelocities() const {
    assert( _hasBeenInitialized );
    BoundaryVector velocity = _regularizer.toBoundary( _baseFlow );
    return velocity;
}

Scalar NonlinearNavierStokes::N(const State& x) const {
    Flux v = CrossProduct( x.q, x.gamma );
    Scalar g = Curl( v );
    return g;
}

Scalar LinearizedNavierStokes::N(const State& x) const {
    Flux v = CrossProduct( _x0.q, x.gamma );
    v += CrossProduct( x.q, _x0.gamma );
    Scalar g = Curl( v );
    return g;
}
    
    
} // namespace ibpm
