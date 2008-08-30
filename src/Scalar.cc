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
// $HeadURL$

#include "Scalar.h"
#include "Flux.h"
#include "Grid.h"
#include <blitz/array.h>
#include <iostream>
using namespace std;

namespace ibpm {

Scalar& Scalar::laplacian(const Scalar& f) {
	// WORK HERE
	return *this;
}

Scalar& Scalar::laplacianInverse(const Scalar& f) {
	// WORK HERE
	return *this;
}

// Print the whole field to standard output
void Scalar::print() const {
    int nx = _grid.getNx();
    int ny = _grid.getNy();
    for(int i = 0; i <= nx; ++i) {
        for (int j=0; j <= ny; ++j) {
            cout << _data(i,j) << " ";
        }
        cout << endl;
    }
    
}

} // namespace
