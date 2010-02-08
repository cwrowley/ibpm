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
    ) :
    _xShift(0.),
	_yShift(0.) {
    resize( nx, ny, ngrid, length, xOffset, yOffset );
}
	
Grid::Grid(
   int nx,
   int ny,
   int ngrid,
   double length,
   double xOffset,
   double yOffset,
   double xShift,
   double yShift
   )
   {
   resize( nx, ny, ngrid, length, xOffset, yOffset );
   setXShift( xShift );
   setYShift( yShift );
}
	
/// Default constructor: set all parameters to zero
Grid::Grid() {
    _nx = 0;
    _ny = 0;
    _ngrid = 0;
    _xOffset = 0;
    _yOffset = 0;
    _dx = 0;
    _xShift = 0.;
	_yShift = 0.;
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
	
/// Set all grid parameters
void Grid::resize(
	int nx,
	int ny,
	int ngrid,
	double length,
	double xOffset,
	double yOffset,
	double xShift,
	double yShift
	) {
	assert( ngrid == 1 || ( nx % 4 == 0 && ny % 4 == 0 ) );
	_nx = nx;
	_ny = ny;
	_ngrid = ngrid;
	_xOffset = xOffset;
	_yOffset = yOffset;
	_xShift = xShift;
	_yShift = yShift;
	_dx = length / nx;
}

// Return the x-coordinate of the left-most gridpoint of level lev
double Grid::getXOffset(int lev) const {
    return _xOffset + 0.5 * ( _xShift - 1 ) * ( ( 1 << lev ) - 1 ) * (_nx * _dx);
}
    
// Return the y-coordinate of the bottom gridpoint of level lev
double Grid::getYOffset(int lev) const {
    return _yOffset + 0.5 * ( _yShift - 1 ) * ( ( 1 << lev ) - 1 ) * (_ny * _dx);
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

// Set the shift parameter in x
void Grid::setXShift( double xShift ) {
    assert( fabs( xShift ) <= 1 );
    assert( fmod( xShift*_nx,4 ) == 0 );
    _xShift = xShift;
}

// Return the shift parameter in x
double Grid::getXShift() const {
    return _xShift;
}

// Set the shift parameter in y
void Grid::setYShift( double yShift ) {
	assert( fabs( yShift ) <= 1 );
	assert( fmod( yShift*_ny,4 ) == 0 );
	_yShift = yShift;
}
	
// Return the shift parameter in y
double Grid::getYShift() const {
	return _yShift;
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

// Compare two grids
bool Grid::isEqualTo( const Grid& grid2 ) const {
    bool nx_eq = ( _nx == grid2.Nx() );
    bool ny_eq = ( _ny == grid2.Ny() );
    bool ngrid_eq = ( _ngrid == grid2.Ngrid() );
    bool dx_eq = ( _dx == grid2.Dx() );
    bool xOffset_eq = ( _xOffset == grid2.getXEdge(0,0) );
    bool yOffset_eq = ( _yOffset == grid2.getYEdge(0,0) );
    bool xShift_eq = ( _xShift == grid2.getXShift() );
    bool yShift_eq = ( _yShift == grid2.getYShift() );
    return( nx_eq * ny_eq * ngrid_eq * dx_eq * xOffset_eq * yOffset_eq * xShift_eq * yShift_eq );
}

} // namespace
