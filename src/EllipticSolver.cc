// EllipticSolver.cc
//
// Description:
// Implementation of the EllipticSolver class
//
// Author(s):
// Clancy Rowley
//
// Date: 30 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "EllipticSolver.h"
#include <math.h>


namespace ibpm {

typedef EllipticSolver2d::Array2d Array2d;

EllipticSolver::EllipticSolver( const Grid& grid ) :
    _ngrid( grid.Ngrid() ),
    _dx( grid.Dx() ),
    _solvers( grid.Ngrid() )
    {}

// Deallocate 2d solvers
EllipticSolver::~EllipticSolver() {
    for (int lev=0; lev<_ngrid; ++lev) {
        delete _solvers[lev];
    }
}
    
    // Create 2d solvers
void EllipticSolver::init() {
    for (int lev=0; lev<_ngrid; ++lev) {
        // calculate grid spacing on this grid level
        double dx = _dx * ( 1 << lev );
        _solvers[lev] = create2dSolver( dx );
    }
}

// Multi-domain elliptic solver
void EllipticSolver::solve( const Scalar& f, Scalar& u ) const {
    assert( f.Ngrid() == _ngrid );
    assert( f.Ngrid() == u.Ngrid() );
    assert( f.Nx() == u.Nx() );
    assert( f.Ny() == u.Ny() );
    
    // First "coarsify" right-hand side (f) to coarse grids
    Scalar rhs = f;
    rhs.coarsify();
    // Solve coarsest grid first, then finer grids
    for (int lev = f.Ngrid() - 1; lev >= 0; --lev ) {
        // Get slices of input and output data at current grid level
        Array2d rhs1 = rhs[lev];
        Array2d u1 = u[lev];
        // if on the coarsest grid, solve with zero bcs
        if (lev == f.Ngrid() - 1) {
            _solvers[lev]->solve( rhs1, u1 );
        }
        else {
            // Get boundary condition from next coarser grid
            BC bc(f.Nx(), f.Ny());
            u.getBC( lev, bc );
            // solve with specified boundary conditions
            _solvers[lev]->solve( rhs1, bc, u1 );
        }
    }
}

Scalar EllipticSolver::solve( const Scalar& f ) const {
    Scalar u( f.getGrid() );
    solve( f, u );
    return u;
}

/******************************************************************************/

PoissonSolver::PoissonSolver( const Grid& grid ) :
    EllipticSolver( grid ),
    _nx( grid.Nx() ),
    _ny( grid.Ny() )
    {
    init();
}

EllipticSolver2d* PoissonSolver::create2dSolver( double dx ) {
    return new PoissonSolver2d( _nx, _ny, dx );
}

/******************************************************************************/

HelmholtzSolver::HelmholtzSolver( const Grid& grid, double alpha ) :
    EllipticSolver( grid ),
    _nx( grid.Nx() ),
    _ny( grid.Ny() ),
    _alpha( alpha )
    {
    init();
}

EllipticSolver2d* HelmholtzSolver::create2dSolver( double dx ) {
    return new HelmholtzSolver2d( _nx, _ny, dx, _alpha );
}
    
} // namespace ibpm
