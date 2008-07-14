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

/// Set *this to the curl of the argument, returning *this.
Flux& Flux::curl(const Scalar& f) {
	assert((f.getGrid()).getNx() == _nx);
  	assert((f.getGrid()).getNy() == _ny);	
 	
	Range all = Range::all();		
	// X direction
	Range J(0, _ny-2);
	_xdata(Range(0,_nx-1), J) = f.getData()(all, J+1) - f.getData()(all,J);
	// Assume zero b.c.s at ghost cells (nx,:), (:, ny).
	_xdata(Range(0, _nx-1), _ny-1) = - f.getData()(all, _ny-1); 
	_xdata(_nx, all) = 0; 
		
	// Y direction
	Range I(0, _nx-2);	
	_ydata(I, Range(0, _ny-1)) = f.getData()(I, all) - f.getData()(I+1, all);
	// Assume zero b.c.s at ghost cells (nx,:), (:, ny).
	_ydata(_nx-1, Range(0, _ny-1)) = f.getData()(_nx-1,all); 
	_ydata(all, _ny) = 0;	
		
	return *this;
}

/// Set *this to the gradient of the argument, returning *this.
Flux& Flux::gradient(const Scalar& f) {
	assert((f.getGrid()).getNx() == _nx);
  	assert((f.getGrid()).getNy() == _ny);
	
	Range all = Range::all();			
	// X direction 	
	Range I(1, _nx-1);
	_xdata(I, all) = f.getData()(I,all) - f.getData()(I-1,all); 
	// Assume zero b.c.s at ghost cells (-1,:), (nx, :).
	_xdata(0, all) = f.getData()(0, all);	
	_xdata(_nx,all) = - f.getData()(_nx-1,all); 	

	// Y direction 	
	Range J(1, _ny-1);
	_ydata(all,J) = f.getData()(all, J) - f.getData()(all, J-1);
	// Assume zero b.c.s at ghost cells (:,-1), (:, ny).
	_ydata(all, 0) = f.getData()(all, 0);		
	_ydata(all, _ny) = - f.getData()(all, _ny-1);	
		
	return *this;
}

/// Return the inner product of *this and the argument
double Flux::dot(const Flux& q) const {
	assert(q._nx == this->_nx);
	assert(q._ny == this->_ny);
	double dp = sum(this->_xdata * q._xdata) + sum(this->_ydata * q._ydata);
	return dp;	
}
