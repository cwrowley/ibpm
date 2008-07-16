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
#include "Flux.h"
#include "Grid.h"
#include <blitz/array.h>

Scalar& Scalar::sinTransform(const Scalar& f) {
	// WORK HERE
	return *this;
}

Scalar& Scalar::sinTransformInv(const Scalar& f) {
	// WORK HERE
	return *this;
}

Scalar& Scalar::laplacian(const Scalar& f) {
	// WORK HERE
	return *this;
}

Scalar& Scalar::laplacianInverse(const Scalar& f) {
	// WORK HERE
	return *this;
}

//dot product of two Scalar products. Only the inner nodes are counted.
double Scalar::dot(const Scalar& f) {
	assert(f._nx == this->_nx);
	assert(f._ny == this->_ny);
	Range I(1, _nx-1);
	Range J(1, _ny-1); 
	double dp = sum(this->_data(I, J) * f._data(I, J));
	dp *=  pow(this->_grid.getDx(),2); 
	return dp;
}
