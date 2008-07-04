// Scalar.cc
//
// Description:
// Implementation of the Scalar class, using Blitz++ arrays
//
// Author(s):
// Clancy Rowley
//
// Date: 4 Jul 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL:// $Header$

#include "Scalar.h"
#include "Grid.h"
#include <blitz/array.h>

Scalar::Scalar(const Grid& grid) {
	_grid = &grid;
	_nx = _grid->getNx();
	_ny = _grid->getNy();
	// allocate memory here
}

Scalar::~Scalar() {
	// deallocate memory here
}

// Scalar f, f1, f2;
// double a;

// TO DO: overload operator f(i,j)

// TODO: overload operator for f1 + f2

// TODO: overload operator for f1 - f2

// TODO: overload operator for a + f

// TODO: overload operator for a - f

// TODO: overload operator for f + a

// TODO: overload operator for f - a

// TODO: overload operator for a * f

// TODO: overload operator for f * a

// TODO: overload operator for f / a

Scalar Scalar::curl(const Flux& q) {
	// WORK HERE
	return *this;
}

Scalar Scalar::sinTransform(const Scalar& f) {
	// WORK HERE
	return *this;
}

Scalar Scalar::sinTransformInv(const Scalar& f) {
	// WORK HERE
	return *this;
}

Scalar Scalar::laplacian(const Scalar& f) {
	// WORK HERE
	return *this;
}

Scalar Scalar::laplacianInverse(const Scalar& f) {
	// WORK HERE
	return *this;
}

Scalar Scalar::divergence(const Flux& q) {
	// WORK HERE
	return *this;
}

double Scalar::dot(const Scalar& f) {
	// WORK HERE
	return 0;
}
