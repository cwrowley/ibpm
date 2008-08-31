// NavierStokesModel.cc
//
// Description:
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
    _grid(grid),
    _geometry(geometry),
    _regularizer( grid, geometry ),
    _linearTermEigenvalues( grid ),
    _eigGammaToStreamfunction( grid ),
    _baseFlow(q_potential),
    _ReynoldsNumber(Reynolds),
    _hasBeenInitialized(false)
{}

void NavierStokesModel::init() {
    if ( _hasBeenInitialized ) return;  // do only once
    // calculate eigenvalues of Laplacian
    int nx = _grid.getNx();
    int ny = _grid.getNy();
    double bydx2 = 1. / ( _grid.getDx() * _grid.getDx() );
    const double pi = 4. * atan(1.);
    Scalar eigLaplacian(_grid);
    eigLaplacian = 1.;
    // Loop over only interior points
    for (int i=1; i < nx; ++i ) {
        for (int j=1; j < ny; ++j ) {
            eigLaplacian(i,j) = 2. * ( cos( (pi * i) / nx ) +
                                       cos( (pi * j) / ny ) - 2. ) * bydx2; 
        }
    }

    // eigenvalues of the linear operation from circulation to
    // streamfunction:
    //    Laplacian psi = - omega
    //    gamma = omega * dx^2
    //    psi = (-Laplacian^{-1} / dx^2) gamma
    _eigGammaToStreamfunction = (-bydx2) / eigLaplacian;
    
    // calculate linear term
    double beta = 1. / _ReynoldsNumber;
    _linearTermEigenvalues = beta * eigLaplacian;

    // Update regularizer
    _regularizer.update();
    
    _hasBeenInitialized = true;
}

NavierStokesModel::~NavierStokesModel() {}

} // namespace ibpm
