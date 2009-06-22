// BC.cc
//
// Description:
// Implementation of the BC class
//
// Author(s):
// Clancy Rowley
//
// Date: 28 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

// Boundary data is stored in a 1d array.
// For an (8 x 4) grid, the data is arranged as follows:
//
//  4  5  6  7  8  9 10 11 12
//  3                      13
//  2                      14
//  1                      15
//  0 23 22 21 20 19 18 17 16
//
// Note that the total number of boundary points is
//   (nx+1) * 2 + (ny-1) * 2 = 2 * (nx + ny)

#include "BC.h"

namespace ibpm {
    
BC::BC( int nx, int ny ) :
    _nx( nx ),
    _ny( ny ),
    _data( 2*(nx+ny) )
    {
    _data = 0.;
}

BC::BC( const BC& bc ) :
    _nx( bc._nx ),
    _ny( bc._ny ),
    _data( 2*(_nx+_ny) ) {

    // copy data
    for (unsigned int i=0; i<_data.Size(); ++i) {
        _data(i) = bc._data(i);
    }
}

BC::~BC()
{}

} // namespace ibpm
