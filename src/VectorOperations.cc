// VectorOperations.cc
//
// Description:
// Definition of VectorOperations functions, using Blitz++ arrays.
//
// Author(s):
// Clancy Rowley, Zhanhua Ma
//
// Date: 15 Jul 2008
//
// $Revision: $
// $LastChangedDate: $
// $LastChangedBy: zma $
// $HeadURL:// $Header$

#include "Scalar.h"
#include "Flux.h"
#include "Grid.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

// Return the curl of Flux q, as a Scalar object
Scalar curl(const Flux& q) { 
	Scalar f(q.getGrid());
	Range I(1, f._nx-1);
	Range J(1, f._ny-1);
	f._data(I, J) = q._ydata(I,J) - q._ydata(I-1,J)  
					+ q._xdata(I,J-1) - q._xdata(I,J);
	//Assume zero b.c.s at the boundary nodes.
	Range all = Range::all(); 
	f._data(0, all) = 0;	
	f._data(all, 0) = 0;
	f._data(f._nx, all) = 0;
	f._data(all,f._ny) = 0;		 				 
	return f;
}

// Return the divergence of Flux q, as a Scalar object
Scalar divergence(const Flux& q) {
//	Scalar f(q.getGrid());
//	int nx = f._nx;
//	int ny = f._ny;	
//	Range I(0, nx-1);
//	Range J(0, ny-1);
//	Range all = Range::all();
//	f._data(all,all) = q._xdata(I+1,J) - q._xdata(I,J)
//							+ q._ydata(I,J+1) - q._ydata(I,J);
//	return f;
}

// Return the curl of Scalar f, as a Flux object.
Flux curl(const Scalar& f) {
	Flux q(f.getGrid());
	Range all = Range::all();		
	// X direction
	Range J(0, q._ny-1);
	q._xdata(all, J) = f._data(all, J+1) - f._data(all,J);			
	// Y direction
	Range I(0, q._nx-1);	
	q._ydata(I, all) = f._data(I, all) - f._data(I+1, all);
			
	return q;
}

// Return the gradient of Scalar f, as a Scalar object.
Flux gradient(const Scalar& f) {
//	Flux q(f.getGrid());
//	int nx = q._nx;
//	int ny = q._ny;	
//	Range all = Range::all();			
//	// X direction 	
//	Range I(1, nx-1);
//	q._xdata(I, all) = f._data(I,all) - f._data(I-1,all); 
//	// Assume zero b.c.s at ghost cells (-1,:), (nx, :).
//	q._xdata(0, all) = f._data(0, all);	
//	q._xdata(nx,all) = - f._data(nx-1,all); 	
//
//	// Y direction 	
//	Range J(1, ny-1);
//	q._ydata(all,J) = f._data(all, J) - f._data(all, J-1);
//	// Assume zero b.c.s at ghost cells (:,-1), (:, ny).
//	q._ydata(all, 0) = f._data(all, 0);		
//	q._ydata(all, ny) = - f._data(all, ny-1);	
//		
//	return q;
}


