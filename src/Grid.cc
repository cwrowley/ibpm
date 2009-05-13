// Grid.cc
//
// Description:
// Implementation of the Grid class
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

#include "Grid.h"
#include <math.h>
#include <assert.h>

namespace ibpm {

Grid::Grid(
    int nx,
    int ny,
    int ngrid,
    double length,
    double xOffset,
    double yOffset
    ) {
    resize( nx, ny, ngrid, length, xOffset, yOffset );
}

/// Default constructor: set all parameters to zero
Grid::Grid() {
    _nx = 0;
    _ny = 0;
    _ngrid = 0;
    _xOffset = 0;
    _yOffset = 0;
    _dx = 0;
};

/// Set all grid parameters
void Grid::resize(
    int nx,
    int ny,
    int ngrid,
    double length,
    double xOffset,
    double yOffset
    ) {
    assert( ngrid == 1 || ( nx % 4 == 0 && ny % 4 == 0 ) );
    _nx = nx;
    _ny = ny;
    _ngrid = ngrid;
    _xOffset = xOffset;
    _yOffset = yOffset;
    _dx = length / nx;
}

// Return the x-coordinate of the left-most gridpoint of level lev
double Grid::getXOffset(int lev) const {
    return _xOffset -  0.5 * ( ( 1 << lev ) - 1 ) * (_nx * _dx);
}
    
// Return the y-coordinate of the bottom gridpoint of level lev
double Grid::getYOffset(int lev) const {
    return _yOffset -  0.5 * ( ( 1 << lev ) - 1 ) * (_ny * _dx);
}
    
    
// Return the x-coordinate of the center of cell i  (i in 0..m-1)
double Grid::getXCenter(int lev, int i) const {
    assert( lev >= 0 && lev <  _ngrid );
    assert(   i >= 0 &&   i <= _nx );
    return getXOffset(lev) + (i+0.5) * Dx(lev);
}

// Return the y-coordinate of the center of cell j  (j in 0..n-1)
double Grid::getYCenter(int lev, int j) const {
    assert( lev >= 0 && lev <  _ngrid );
    assert(   j >= 0 &&   j <= _ny );
    return getYOffset(lev) + (j+0.5) * Dx(lev);
}

// Return the x-coordinate of the left edge of cell i  (i in 0..m)
double Grid::getXEdge(int lev, int i) const {
    assert( lev >= 0 && lev <  _ngrid );
    assert(   i >= 0 &&   i <= _nx );
    return getXOffset(lev) + i * Dx(lev);
}

// Return the y-coordinate of the bottom edge of cell j  (j in 0..n)
double Grid::getYEdge(int lev, int j) const {
    assert( lev >= 0 && lev <  _ngrid );
    assert(   j >= 0 &&   j <= _ny );
    return getYOffset(lev) + j * Dx(lev);
}

// Return the grid index i corresponding to the given x-coordinate 
// Currently, only works for the finest grid level.
int Grid::getXGridIndex( double x ) const {
	double xpos = x - _xOffset;
	assert ( xpos <=  _dx * _nx);
	int i = int(floor( xpos / _dx ));
	if ( (xpos - i * _dx) >= (_dx / 2) ) {
		 i = int(ceil( xpos / _dx ));
	}
	return i;
}

// Return the grid index j corresponding to the given j-coordinate 
// Currently, only works for the finest grid level.
int Grid::getYGridIndex( double y ) const {
	double ypos = y - _yOffset;
	assert ( ypos <=  _dx * _ny);
	int j = int(floor( ypos / _dx ));
	if ( (ypos - j * _dx) >= (_dx / 2) ) {
		 j = int(ceil( ypos / _dx ));
	}
	return j;
}

} // namespace

