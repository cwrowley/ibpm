#ifndef _SCALAR_H_
#define _SCALAR_H_

#include "Array.h"
#include "Field.h"
#include "Grid.h"
#include "Flux.h"

namespace ibpm {

/*!
    \file Scalar.h
    \class Scalar

    \brief Store a 2D array of scalar values, located at cell nodes.
    
    For a grid with nx cells in the x-direction and ny cells in the y-direction,
    there are (nx+1)*(ny+1) nodes.  
    
    There are (nx-1)*(ny-1) inner nodes, and 2*(nx+ny) boundary nodes.
    
    Also provides scalar-valued math operations, such as discrete sin transform, 
    and Laplacian and its inverse, of scalar fields.

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

    // Print the whole field to standard output
    void print() const;
    
    /// Copy assignment
    inline Scalar& operator=(const Scalar& f) {
        assert( f.Nx() == Nx() );
        assert( f.Ny() == Ny() );
        _data = f._data;
        return *this;
    }

    /// Copy assignment from double
    inline Scalar& operator=(double a) {
        _data = a;
        return *this;
    }

    /// f(i,j) refers to the value at index (i,j)
    inline double& operator()(int i, int j) {
        assert(i>=0  && i<= Nx() );
        assert(j>=0  && j<= Ny() );
        return _data(i,j);
    }
    
    /// f(i,j) refers to the value at index (i,j)
    inline double operator()(int i, int j) const{
        assert(i>=0  && i<= Nx() );
        assert(j>=0  && j<= Ny() );
        return _data(i,j);
    }
    
    /// f += g
    inline Scalar& operator+=(const Scalar& f) {
        assert(f.Nx() == Nx() );
        assert(f.Ny() == Ny() );
        _data += f._data;
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
        assert( f.Nx() == Nx() );
        assert( f.Ny() == Ny() );
        _data -= f._data;
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
    inline Scalar operator-(Scalar& f) {
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
    inline Scalar operator/(Scalar& f) {
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
    Array::Array2<double> _data;
};

/// -f
inline Scalar operator-(Scalar& f) {
    Scalar g = f;
    g *= -1;
    return g;
}

/// a + f
inline Scalar operator+(double a, Scalar& f) {
    Scalar g = f;
    g += a;
    return g;
}
    
/// a - f
inline Scalar operator-(double a, Scalar& f) {
    Scalar g = -f;
    g += a;
    return g;
}
    
/// a * f
inline Scalar operator*(double a, Scalar& f) {
    Scalar g = f;
    g *= a;
    return g;
}

/// a / f
inline Scalar operator/(double a, Scalar& f) {
    Scalar g(f.getGrid());
    g = a;
    g /= f;
    return g;
}

} // namespace ibpm

#endif /* _SCALAR_H_ */
