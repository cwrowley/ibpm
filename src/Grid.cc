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
    _nx = nx;
    _ny = ny;
    _ngrid = ngrid;
    _xOffset = xOffset;
    _yOffset = yOffset;
    _dx = length / nx;
}

// Return the x-coordinate of the center of cell i  (i in 0..m-1)
double Grid::getXCenter(int i) const {
    assert(i >= 0);
    assert(i < _nx);
    return _xOffset + (i+0.5) * _dx;
}

// Return the y-coordinate of the center of cell j  (j in 0..n-1)
double Grid::getYCenter(int j) const {
    assert(j >= 0);
    assert(j < _ny);
    return _yOffset + (j+0.5) * _dx;
}

// Return the x-coordinate of the left edge of cell i  (i in 0..m)
double Grid::getXEdge(int i) const {
    assert(i >= 0);
    assert(i <= _nx);
    return _xOffset + i * _dx;
}

// Return the y-coordinate of the bottom edge of cell j  (j in 0..n)
double Grid::getYEdge(int j) const {
    assert(j >= 0);
    assert(j <= _ny);
    return _yOffset + j * _dx;
}

// Return the area of the domain
double Grid::getArea() const {
    return _dx * _dx * _nx * _ny;
}

} // namespace

