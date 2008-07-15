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

Scalar& Scalar::curl(const Flux& q) {
	assert((q.getGrid()).getNx() == _nx);
  	assert((q.getGrid()).getNy() == _ny);
	Range I(1, _nx-1);
	Range J(1, _ny-1);
	_data(I, J) = q.getDataY()(I,J) - q.getDataY()(I-1,J)  
					+ q.getDataX()(I,J-1) - q.getDataX()(I,J);
	//Assume zero b.c.s at bottom and left "ghost edges".
	Range all = Range::all(); 
	_data(0, J) = q.getDataY()(0, J) + q.getDataX()(0, J-1) - q.getDataX()(0,J);	
	_data(I, 0) = q.getDataY()(I, 0) - q.getDataY()(I-1, 0) - q.getDataX()(I,0);
	_data(0, 0) = q.getDataY()(0, 0) - q.getDataX()(0,0);				 				 
	return *this;
}

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

Scalar& Scalar::divergence(const Flux& q) {
	assert((q.getGrid()).getNx() == _nx);
  	assert((q.getGrid()).getNy() == _ny);
	Range I(0, _nx-1);
	Range J(0, _ny-1);
	Range all = Range::all();
	_data(all,all) = q.getDataX()(I+1,J) - q.getDataX()(I,J)
							+ q.getDataY()(I,J+1) - q.getDataY()(I,J);	

	return *this;
}

double Scalar::dot(const Scalar& f) {
	assert(f._nx == this->_nx);
	assert(f._ny == this->_ny);
	double dp = sum(this->_data * f._data);
	dp *=  pow(this->_grid.getDx(),2); 
	return dp;
}
