#ifndef _SCALAR_H_
#define _SCALAR_H_

#include "Array.h"
#include "BC.h"
#include "Field.h"
#include "Grid.h"
#include <iostream>  
 
namespace ibpm {

/*!
    \file Scalar.h
    \class Scalar

    \brief Store a 2D array of scalar values, located at interior cell nodes.
    
    For a grid with nx cells in the x-direction and ny cells in the y-direction,
    there are (nx+1)*(ny+1) nodes.  
    
    There are (nx-1)*(ny-1) inner nodes, and 2*(nx+ny) boundary nodes.
    Only the interior nodes are stored in a Scalar, and the boundary nodes are
    always zero.
    
    \author Clancy Rowley
    \author $LastChangedBy$
    \date  4 Jul 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class Scalar : public Field {
public:
    /// Allocate memory for the 2D array
    Scalar( const Grid& grid );

    /// Default constructor: do not allocate memory yet
    Scalar();
    
    /// Allocate new array, copy the data
    Scalar( const Scalar& f );
    
    /// Destructor
    ~Scalar();

    /// Reassign the grid parameters and allocate memory based on the new grid
    void resize( const Grid& grid );

    /// Print the whole field to standard output
    void print() const;
    
    /// "Coarsify" the Scalar quantity
    ///  - Fine grid is left unchanged
    ///  - Coarse values that correspond to points on the fine grid are replaced by
    ///    averaged value of fine gridpoints
    void coarsify();
    
    /// Copy assignment
    inline Scalar& operator=(const Scalar& f) {
        assert( f.Ngrid() == Ngrid() );
		assert( f.Nx() == Nx() );
        assert( f.Ny() == Ny() );
        for (unsigned int i=0; i<_data.Size(); ++i) {
            _data(i) = f._data(i);
        }
        return *this;
    }

    /// Copy assignment from double
    inline Scalar& operator=(double a) {
        _data = a;
        return *this;
    }

    /// f(i,j) refers to the value at index (i,j)
    inline double& operator()(int lev, int i, int j) {
        assert( lev >= 0 && lev < Ngrid() );
        assert( i>=1  && i< Nx() );
        assert( j>=1  && j< Ny() );
        return _data(lev,i,j);
    }
    
    /// f(i,j) refers to the value at index (i,j)
    inline double operator()(int lev, int i, int j) const{
        assert( lev >= 0 && lev < Ngrid() );
        assert( i>=1  && i< Nx() );
        assert( j>=1  && j< Ny() );
        return _data(lev,i,j);
    }
		
    /// f[lev] returns a 2d array of grid level lev
    inline Array::Array2<double> operator[](int lev) {
        return _data[lev];
    }
        
    /// f[lev] returns a 2d array of grid level lev
    inline const Array::Array2<double> operator[](int lev) const {
        return _data[lev];
    }
    
    /// \brief Compute the boundary values at level lev from the next coarser
    /// grid.
    /// \param[in] lev is the grid level for which the bounday values are
    ///             desired; must be in the range 0..Ngrid-2
    /// \param[out] bc contains the boundary values computed
    void getBC( int lev, BC& bc ) const;
    
    /// f += g
    inline Scalar& operator+=(const Scalar& f) {
        assert( f.Ngrid() == Ngrid() );
        assert(f.Nx() == Nx() );
        assert(f.Ny() == Ny() );
        for ( unsigned int i=0; i<_data.Size(); ++i ) {
            _data(i) += f._data(i);
        }
        return *this;
    }

    /// f += a
    inline Scalar& operator+=(double a) {
        for ( unsigned int i=0; i<_data.Size(); ++i ) {
            _data(i) += a;
        }
        return *this;
    }

    /// f -= g
    inline Scalar& operator-=(const Scalar& f) {
        assert( f.Ngrid() == Ngrid() );
        assert(f.Nx() == Nx() );
        assert(f.Ny() == Ny() );
        for ( unsigned int i=0; i<_data.Size(); ++i ) {
            _data(i) -= f._data(i);
        }
        return *this;
    }

    /// f -= a
    inline Scalar& operator-=(double a) {
        for ( unsigned int i=0; i<_data.Size(); ++i ) {
            _data(i) -= a;
        }
        return *this;
    }

    /// f + g
    inline Scalar operator+(const Scalar& f) {
        assert( f.Ngrid() == Ngrid() );
        assert( f.Nx() == Nx() );
        assert( f.Ny() == Ny() );
        Scalar g = *this;
        g += f;
        return g;
    }
    
    /// f + a
    inline Scalar operator+(double a) {
        Scalar g = *this;
        g += a;
        return g;
    }
    
    /// f - g
    inline Scalar operator-(const Scalar& f) {
        assert( f.Ngrid() == Ngrid() );
        assert( f.Nx() == Nx() );
        assert( f.Ny() == Ny() );
        Scalar g = *this;
        g -= f;
        return g;
    }
    
    /// f - a
    inline Scalar operator-(double a) {
        Scalar g = *this;
        g -= a;
        return g;
    }
    
    /// f *= g
    inline Scalar& operator*=(const Scalar& f) {
        assert( f.Ngrid() == Ngrid() );
        assert( f.Nx() == Nx() );
        assert( f.Ny() == Ny() );
        for ( unsigned int i=0; i<_data.Size(); ++i) {
            _data(i) *= f._data(i);
        }
        return *this;
    }

    /// f *= a
    inline Scalar& operator*=(double a) {
        _data *= a;
        return *this;
    }

    /// f * g
    inline Scalar operator*(const Scalar& f) {
        assert( f.Ngrid() == Ngrid() );
        assert( f.Nx() == Nx() );
        assert( f.Ny() == Ny() );
        Scalar g = *this;
        g *= f;
        return g;
    }
    
    /// f * a
    inline Scalar operator*(double a) {
        Scalar g = *this;
        g *= a;
        return g;
    }
    
    /// f /= g
    inline Scalar& operator/=(const Scalar& f) {
        assert( f.Ngrid() == Ngrid() );
        assert( f.Nx() == Nx() );
        assert( f.Ny() == Ny() );
        _data /= f._data;
        return *this;
    }

    /// f /= a
    inline Scalar& operator/=(double a) {
        _data /= a;
        return *this;
    }

    /// f / g
    inline Scalar operator/(const Scalar& f) {
        assert( f.Ngrid() == Ngrid() );
        assert( f.Nx() == Nx() );
        assert( f.Ny() == Ny() );
        Scalar g = *this;
        g /= f;
        return g;
    }
    
    /// f / a
    inline Scalar operator/(double a) {
        Scalar g = *this;
        g /= a;
        return g;
    }


private:
    Array::Array3<double> _data;
};

/// -f
inline Scalar operator-(const Scalar& f) {
    Scalar g = f;
    g *= -1;
    return g;
}

/// a + f
inline Scalar operator+(double a, const Scalar& f) {
    Scalar g = f;
    g += a;
    return g;
}
    
/// a - f
inline Scalar operator-(double a, const Scalar& f) {
    Scalar g = -f;
    g += a;
    return g;
}
    
/// a * f
inline Scalar operator*(double a, const Scalar& f) {
    Scalar g = f;
    g *= a;
    return g;
}

/// a / f
inline Scalar operator/(double a, const Scalar& f) {
    Scalar g(f.getGrid());
    g = a;
    g /= f;
    return g;
}

} // namespace ibpm

#endif /* _SCALAR_H_ */
