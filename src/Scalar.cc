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

Scalar::Scalar( const Grid& grid ){
    resize( grid );
}

/// Default constructor: do not allocate memory yet
Scalar::Scalar() {
    _nx = 0;
    _ny = 0;
}
    
/// Allocate new array, copy the data
Scalar::Scalar( const Scalar& f ) {
    resize( f._grid );
    // copy data
    _data = f._data;
}
    
/// Deallocate memory in the destructor
Scalar::~Scalar() {
    // deallocation automatic for Blitz++ arrays
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

void Scalar::resize( const Grid& grid ) {
    _grid = grid;
    _nx = grid.getNx();
    _ny = grid.getNy();
    _data.resize( _nx+1, _ny+1 );
}


} // namespace
