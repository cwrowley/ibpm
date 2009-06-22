#ifndef _BC_H_
#define _BC_H_

#include "Array.h"

namespace ibpm {
    
/*!
    \file BC.h
    \class BC

    \brief Class for storing and accessing boundary conditions for a 2d
    scalar field.

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 28 Sep 2008
    \date $LastChangedDate$
    \version $Revision$
*/

// For an (8 x 4) grid, the data is arranged as follows:
//
//  4  5  6  7  8  9 10 11 12
//  3                      13
//  2                      14
//  1                      15
//  0 23 22 21 20 19 18 17 16

class BC {
public:
    /// \brief Constructor: initializes boundary data to zero
    /// \param[in] nx number of cells in x-direction
    /// \param[in] ny number of cells in y-direction
    BC( int nx, int ny );
    BC( const BC& bc );
    ~BC();
    
    /// \brief Return reference to value on left boundary, at index j
    /// \param[in] j must be from 0..ny
    inline double& left( int j ) {
        assert( j >= 0 && j <= _ny );
        return _data(j);
    }
    inline const double& left( int j ) const {
        assert( j >= 0 && j <= _ny );
        return _data(j);
    }
    
    /// \brief Return reference to value on right boundary, at index j
    /// \param[in] j must be from 0..ny
    inline double& right( int j ) {
        assert( j >= 0 && j <= _ny );
        return _data( 2 * _ny + _nx - j);
    }
    inline const double& right( int j ) const {
        assert( j >= 0 && j <= _ny );
        return _data( 2 * _ny + _nx - j);
    }
    
    /// \brief Return reference to value on top boundary, at index i
    /// \param[in] i must be from 0..nx
    inline double& top( int i ) {
        assert( i >= 0 && i <= _nx );
        return _data( _ny + i );
    }
    inline const double& top( int i ) const {
        assert( i >= 0 && i <= _nx );
        return _data( _ny + i );
    }
    
    /// \brief Return reference to value on bottom boundary, at index i
    /// \param[in] j must be from 0..nx
    inline double& bottom( int i ) {
        assert( i >= 0 && i <= _nx );
        if (i == 0) return _data(0);
        else return _data( 2 * (_nx + _ny) - i );
    }
    inline const double& bottom( int i ) const {
        assert( i >= 0 && i <= _nx );
        if (i == 0) return _data(0);
        else return _data( 2 * (_nx + _ny) - i );
    }
    
    inline BC& operator=(double a) {
        _data = a;
        return *this;
    }
    
    inline int Nx() const { return _nx; }
    inline int Ny() const { return _ny; }
    
private:
    int _nx;
    int _ny;
    Array::Array1<double> _data;
};

} // namespace ibpm

#endif /* _BC_H_ */
