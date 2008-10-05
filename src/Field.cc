// Field.cc
//
// Description:
// Implementation of the Field class (base class for Scalar, Flux)
//
// Author(s):
// Clancy Rowley
//
// Date: 6 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "Grid.h"
#include "Field.h"

namespace ibpm {

Field::Field() {
    // Note: cannot set nx to zero or computation of dx will divide by zero
    // Set to -1 to indicate no grid defined
    int nx = -1;
    int ny = -1;
    int ngrid = 1;
    double length = 0;
    double xOffset = 0;
    double yOffset = 0;
    _grid.resize( nx, ny, ngrid, length, xOffset, yOffset );
}

Field::Field(const Grid& grid) {
    _grid = grid;
}

Field::Field(const Field& field) {
    _grid = field._grid;
}

Field::~Field() {}

} // namespace ibpm
