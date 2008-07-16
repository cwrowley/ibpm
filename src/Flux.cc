// Flux.cc
//
// Description:
// Implementation of the Flux class, using Blitz++ arrays
//
// Author(s):
// Clancy Rowley
//
// Date: 4 Jul 2008
//
// $Revision: 24 $
// $LastChangedDate: 2008-07-04 00:50:47 -0400 (Fri, 04 Jul 2008) $
// $LastChangedBy: clancy $
// $HeadURL:// $Header$

#include "Flux.h"
#include "Scalar.h"
#include "Grid.h"
#include <blitz/array.h>

// dot product of *this and the argument. Only the "inner fluxes" are counted.
double Flux::dot(const Flux& q) const {
	assert(q._nx == this->_nx);
	assert(q._ny == this->_ny);
	Range all = Range::all();
	Range I(1,_nx-1);
	Range J(1,_ny-1);
	double dp = sum(this->_xdata(I,J) * q._xdata(I,J)); 
	dp +=  sum(this->_ydata(I,J) * q._ydata(I,J));
	dp *=  pow(this->_grid.getDx(),2); 
	return dp;	
}
