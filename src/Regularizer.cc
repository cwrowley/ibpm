// Scalar.cc
//
// Description:
// Implementation of the Regularizer class
//
// Author(s):
// Clancy Rowley
//
// Date: 7 Jul 2008
//
// $Revision: 28 $
// $LastChangedDate: 2008-07-04 01:35:13 -0400 (Fri, 04 Jul 2008) $
// $LastChangedBy: clancy $
// $HeadURL:// $Header$

#include "Regularizer.h"

Regularizer(const Grid& grid, const Geometry& geometry) : 
	_grid(grid),
	_geometry(geom)
{
	// TODO: allocate memory for internal data structures
}

void update() {
	// TODO: work here
}

Scalar toGrid(const BoundaryVector& u1) {
	Scalar u2(_grid);
	u2 = 0;
	// TODO: work here
	return u2;
}

BoundaryVector toBoundary(const Scalar& u2) {
	BoundaryVector u1(_geometry.getNumPoints());
	u1 = 0;
	// TODO: work here
	return u1;
}

