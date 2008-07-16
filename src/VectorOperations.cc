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

	// Compute curl at interior nodes
	Range I(1, f._nx-1);
	Range J(1, f._ny-1);
	f._data(I, J) = q._ydata(I,J) - q._ydata(I-1,J)  
					+ q._xdata(I,J-1) - q._xdata(I,J);

	// Return zero values at the boundary nodes.
	Range all = Range::all(); 
	f._data(0, all) = 0;	
	f._data(all, 0) = 0;
	f._data(f._nx, all) = 0;
	f._data(all,f._ny) = 0;		 				 
	return f;
}

// Return the curl of Scalar f, as a Flux object.
Flux curl(const Scalar& f) {
	Flux q(f.getGrid());
	Range all = Range::all();		

	// Compute X-component of curl
	Range J(0, q._ny-1);
	q._xdata(all, J) = f._data(all, J+1) - f._data(all,J);			

	// Compute Y-component of curl
	Range I(0, q._nx-1);	
	q._ydata(I, all) = f._data(I, all) - f._data(I+1, all);
			
	return q;
}


