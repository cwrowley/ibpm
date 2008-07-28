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
    double dp = 0;
    // Add x-fluxes
    for (int i=1; i<q._nx; ++i) {
        for (int j=1; j<q._ny; ++j) {
            dp += q(X,i,j) * (*this)(X,i,j);
        }
    }
    
    // Add y-fluxes
    for (int i=1; i<q._nx; ++i) {
        for (int j=1; j<q._ny; ++j) {
            dp += q(Y,i,j) * (*this)(Y,i,j);
        }
    }
    double dx = this->_grid.getDx();
    return dp * dx * dx;
}
