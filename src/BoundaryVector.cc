// BoundaryVector.cc
//
// Description:
// Implementation of the BoundaryVector class, using Blitz++ arrays
//
// Author(s):
// Clancy Rowley
//
// Date: 7 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "BoundaryVector.h"

namespace ibpm {

BoundaryVector::BoundaryVector() {}

BoundaryVector::BoundaryVector(int numPoints ) {
    resize( numPoints );
}

/// Allocate a new BoundaryVector, copy the data
BoundaryVector::BoundaryVector(const BoundaryVector& f) {
    resize( f._numPoints );
    _data = f._data;
}

void BoundaryVector::resize( int numPoints ) {
    _numPoints = numPoints;
    // blitz: _data.resize( _numPoints * XY );
    _data.Deallocate();
    _data.Allocate( _numPoints * XY );
}

} // namespace ibpm
