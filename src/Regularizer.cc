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
#include <math>

Regularizer(const Grid& grid, const Geometry& geometry) : 
	_grid(grid),
	_geometry(geom)
{
	// TODO: allocate memory for internal data structures
}

// Update list of relationships between boundary points and cells, and the
// corresponding weights
void Regularizer::update() {
    // TODO: work here
    // For each direction (x and y)
    //     For each point on the boundary
    //         Clear the list of associated cells
    //         For each cell
    //             Find distance between boundary point and cell
    //             If cell is within the radius of support of the delta function
    //                 Add to list of associated cells
    //                 Compute the weight factor

    // Note: use function templates to handle x and y directions?
    // First code in one direction without templates.
}

Flux Regularizer::toGrid(const BoundaryVector& u1) {
	Flux u2(_grid);
	u2 = 0;
    // TODO: work here
    // Allocate a new Flux field, initialized to zero
    // For each direction (x and y)
    //     For each point on the boundary
    //         For each associated cell
    //             Add the weight factor times the boundary value to the flux field
    // Return the new flux field
    return u2;
}

BoundaryVector Regularizer::toBoundary(const Flux& u2) {
	BoundaryVector u1(_geometry.getNumPoints());
	u1 = 0;
    // TODO: work here
    // Allocate a new BoundaryVector, initialized to zero
    // For each direction (x and y)
    //     For each point on the boundary
    //         For each associated cell
    //             Add the weight factor times the flux value to the boundary value
    // Return the new BoundaryVector
	return u1;
}

// Return the value of the regularized delta function
// From Roma, Peskin, and Berger, JCP 1999, eq (22)
double Regularizer::deltaFunction(double x) {
	double h = _grid.getDx();
	double r = abs(x) / h;
	if (r > 1.5) {
		return 0;
	}
	else {
		if (r <= 0.5)
			return (1 + sqrt(1 - 3*r*r))/(3*h);
		else
			return (5 - 3*r - sqrt(1 - 3*(1-r)*(1-r)) ) / (6*h);
	}
}