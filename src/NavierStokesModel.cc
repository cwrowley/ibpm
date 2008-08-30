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
// $HeadURL:// $Header$

#include "NavierStokesModel.h"


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
    _ReynoldsNumber(Reynolds)
    {
    
    // calculate eigenvalues of Laplacian
    int nx = grid.getNx();
    int ny = grid.getNy();    
    const double pi = 4. * atan(1.);
    Scalar eigLaplacian(grid);
    eigLaplacian = 1.;
    // Loop over only interior points
    for (int i=1; i < nx; ++i ) {
        for (int j=1; j < ny; ++j ) {
            eigLaplacian(i,j) = 2. * ( cos( (pi * i) / nx ) +
                                       cos( (pi * j) / ny ) - 2. ) / 
                                ( _grid.getDx() * _grid.getDx() );
        }
    }

    // eigenvalues of the linear operation from circulation to
    // streamfunction:
    //    Laplacian psi = - omega
    //    gamma = omega * dx^2
    //    psi = (-Laplacian^{-1} / dx^2) gamma
    _eigGammaToStreamfunction = -1/( _grid.getDx() * _grid.getDx() ) /
                                eigLaplacian;
    
    // calculate linear term
    double beta = 1. / Reynolds;
    _linearTermEigenvalues = beta * eigLaplacian;

    // Update regularizer
    _regularizer.update();
}

NavierStokesModel::~NavierStokesModel() {}
