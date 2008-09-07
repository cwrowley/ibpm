// Flux.cc
//
// Description:
// Implementation of the Flux class, using Blitz++ arrays
//
// Author(s):
// Clancy Rowley
//
// Date: 4 Jul 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "Flux.h"
#include "Scalar.h"
#include "Grid.h"
#include <blitz/array.h>

namespace ibpm {

Flux::Flux() {}

Flux::Flux( const Grid& grid ) {
    resize( grid );
}

/// Allocate new array, copy the data
Flux::Flux( const Flux& q ) {
    resize( q._grid );
    _data = q._data;
}

/// Set all parameters and reallocate arrays based on the Grid dimensions
void Flux::resize( const Grid& grid ) {
    _grid = grid;
    _nx = grid.getNx();
    _ny = grid.getNy();
    _numXFluxes = _nx * _ny + _ny;
    _numFluxes = 2 * _nx * _ny + _nx + _ny;
    _data.resize( _numFluxes );
}

Flux::~Flux() {} // deallocation automatic for Blitz++ arrays

// Print the X and Y components to standard out (for debugging)
void Flux::print() {
    cout << "X:" << endl;
    for (int i=0; i<=_nx; ++i) {
        for (int j=0; j<_ny; ++j) {
            cout << (*this)(X,i,j) << " ";
        }
        cout << endl;
    }
    cout << "Y:" << endl;
    for (int i=0; i<_nx; ++i) {
        for (int j=0; j<=_ny; ++j) {
            cout << (*this)(Y,i,j) << " ";
        }
        cout << endl;
    }
}

/// Return Flux for a uniform flow with the specified magnitude and dir
Flux Flux::UniformFlow(
    const Grid& grid,
    double magnitude,
    double angle
    ) {
    double dx = grid.getDx();
    double u = magnitude * cos( angle ) * dx;
    double v = magnitude * sin( angle ) * dx;
    
    Flux q(grid);
    Flux::index ind;
    for (ind = q.begin(X); ind != q.end(X); ++ind) {
        q(ind) = u;
    }
    for (ind = q.begin(Y); ind != q.end(Y); ++ind) {
        q(ind) = v;
    }
    return q;
}

} // namespace ibpm
