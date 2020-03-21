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

using std::cout;

namespace ibpm {

Flux::Flux() {}

Flux::Flux( const Grid& grid ) :
    Field( grid ) {
    resize( grid );
}

/// Allocate new array, copy the data
Flux::Flux( const Flux& q ) :
    Field( q ) {
    resize( q.getGrid() );
    // copy the data
    for (unsigned int i=0; i<_data.Size(); ++i) {
        _data(i) = q._data(i);
    }
}

/// Set all parameters and reallocate arrays based on the Grid dimensions
void Flux::resize( const Grid& grid ) {
    setGrid( grid );
    int nx = Nx();
    int ny = Ny();
    _numXFluxes = nx * ny + ny;
    _numFluxes = 2 * nx * ny + nx + ny;
    _data.Deallocate();
    _data.Allocate( Ngrid(), _numFluxes );
}

Flux::~Flux() {} // deallocation automatic for Blitz++ arrays

// Print the X and Y components to standard out (for debugging)
void Flux::print() {
    cout << "X:" << endl;
    for (int lev=0; lev<Ngrid(); ++lev) {
        for (int j=Ny()-1; j>=0; --j) {
            for (int i=0; i<=Nx(); ++i) {
                cout << (*this)(lev,X,i,j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << "Y:" << endl;
    for (int lev=0; lev<Ngrid(); ++lev) {
        for (int j=Ny(); j>=0; --j) {
        for (int i=0; i<Nx(); ++i) {
                cout << (*this)(lev,Y,i,j) << " ";
            }
            cout << endl;
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
    double u = magnitude * cos( angle );
    double v = magnitude * sin( angle );
    
    Flux q(grid);
    Flux::index ind;
    for (int lev=0; lev<grid.Ngrid(); ++lev) {
        double dx = grid.Dx(lev);
        for (ind = q.begin(X); ind != q.end(X); ++ind) {
            q(lev,ind) = u * dx;
        }
        for (ind = q.begin(Y); ind != q.end(Y); ++ind) {
            q(lev,ind) = v * dx;
        }
    }
    return q;
}

void Flux::setFlow(
    TangentSE2 g,
    double xCenter,
    double yCenter
    ) {
    /*
        This function computes the unsteady base flow motion 
          corresponding to a moving body with motion g\in TSE(2)
        The components xdot, ydot, thetadot are calculated in BaseFlow::moveFlow
        The calculation is similar to Flux::UniformFlow, but with a rotational component
    */
    double xdot, ydot, thetadot;
    double xdiff, ydiff;
    g.getVelocity(xdot,ydot,thetadot);
    Flux::index ind;
    for (int lev=0; lev<Ngrid(); ++lev) {
        double dx = Dx(lev);
        for(ind = begin(X); ind != end(X); ++ind) {
            xdiff = x(lev,ind) - xCenter;
            ydiff = y(lev,ind) - yCenter;
            _data(lev,ind) = (xdot -thetadot*ydiff)*dx;
        }
        for(ind = begin(Y); ind != end(Y); ++ind) {
            xdiff = x(lev,ind) - xCenter;
            ydiff = y(lev,ind) - yCenter;
            _data(lev,ind) = (ydot + thetadot*xdiff)*dx;
        }
    }
}

} // namespace ibpm
