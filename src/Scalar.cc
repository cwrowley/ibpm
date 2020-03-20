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
#include <iostream>

using std::cout;

namespace ibpm {

Scalar::Scalar( const Grid& grid ) :
    Field( grid ) {
    resize( grid );
}

/// Default constructor: do not allocate memory yet
Scalar::Scalar() {}
    
/// Allocate new array, copy the data
Scalar::Scalar( const Scalar& f ) :
    Field( f ) {
    resize( f.getGrid() );
    // copy data
    for (unsigned int i=0; i<_data.Size(); ++i) {
        _data(i) = f._data(i);
    }
}
    
/// Deallocate memory in the destructor
Scalar::~Scalar() {
    // deallocation automatic for Arrays
}
    
void Scalar::coarsify() {
    // Fine grid unchanged: start with next finest grid
    for (int lev=1; lev<Ngrid(); ++lev) {
        // Loop over interior gridpoints, that correspond to finer grid
        for (int i=NxExt()+1; i<Nx()/2+NxExt(); ++i) {
            for (int j=NyExt()+1; j<Ny()/2+NyExt(); ++j) {
                // get corresponding point on fine grid
				int ii,jj;
                getGrid().c2f(i,j,ii,jj);
                _data(lev,i,j) = 0.25 * _data(lev-1,ii,jj) +
                    0.125 * ( _data(lev-1,ii+1,jj) + _data(lev-1,ii,jj+1) +
                              _data(lev-1,ii-1,jj) + _data(lev-1,ii,jj-1) ) +
                    0.0625 * ( _data(lev-1,ii+1,jj+1) + _data(lev-1,ii+1,jj-1) +
                               _data(lev-1,ii-1,jj+1) + _data(lev-1,ii-1,jj-1) );
            }
        }
    }
}

// Print the whole field to standard output
void Scalar::print() const {
    int nx = Nx();
    int ny = Ny();
    for( int lev=0; lev<Ngrid(); ++lev ) {
        for (int j=ny-1; j > 0; --j) {
            for(int i = 1; i < nx; ++i) {
                cout << _data(lev,i,j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    
}

void Scalar::resize( const Grid& grid ) {
    setGrid( grid );
    _data.Deallocate();
    // Allocate arrays for interior points:
    //    lev in 0..lev-1
    //    i   in 1..nx-1
    //    j   in 1..ny-1
    _data.Allocate( Ngrid(), Nx() - 1, Ny() - 1, 0, 1, 1 );
}
    
void Scalar::getBC( int lev, BC& bc ) const {
    assert( Nx() == bc.Nx() );
    assert( Ny() == bc.Ny() );
    assert( lev >= 0 && lev < Ngrid()-1 );

    // top and bottom boundaries
	/* if grid is shifted completely up or down, 
	   then all poinsts on the shared boundary must take a value of 0,
	   as required on the boundary of the outermost grid */ 
    for (int i=0; i<=Nx(); ++i) {
        int ii,jj;  // indices on coarse grid
        getGrid().f2c(i,0,ii,jj);
        // copy point that coincides with coarse point
		bc.bottom(i) = (*this)( lev+1, ii, jj );
		bc.top(i) = (*this)( lev+1, ii, Ny()/2+jj );
		
        if ( ++i <= Nx() ) {
            // interpolate point in between coarse points
			bc.bottom(i) = 0.5 * ( (*this)( lev+1, ii, jj ) + (*this)( lev+1, ii+1, jj) );
			bc.top(i) = 0.5 * ( (*this)( lev+1, ii, Ny()/2+jj ) + (*this)( lev+1, ii+1, Ny()/2+jj) );
        }
    }

    // left and right boundaries
    for (int j=0; j<=Ny(); ++j) {
        int ii,jj;  // indices on coarse grid
        getGrid().f2c(0,j,ii,jj);
        // copy point that coincides with coarse point
		bc.left(j) = (*this)( lev+1, ii, jj );
		bc.right(j) = (*this)( lev+1, Nx()/2+ii, jj );
        if ( ++j <= Ny() ) {
            // interpolate point in between coarse points
			bc.left(j) = 0.5 * ( (*this)( lev+1, ii, jj ) + (*this)( lev+1, ii, jj+1) );
			bc.right(j) = 0.5 * ( (*this)( lev+1, Nx()/2+ii, jj ) + (*this)( lev+1, Nx()/2+ii, jj+1 ) );
        }
    }
}



} // namespace
