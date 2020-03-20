#ifndef _FLUX_H_
#define _FLUX_H_

#include "Grid.h"
#include "Field.h"
#include "Direction.h"
#include "Array.h"
#include "TangentSE2.h"
#include <math.h>
#include <iostream>

namespace ibpm {

/*!
    \file Flux.h
    \class Flux
    
    \brief Store a 2D array of fluxes, located at cell edges
    
    For a grid with nx cells in the x-direction and ny cells in the 
    y-direction, there are (nx+1,ny) fluxes in the x-direction, and (nx,ny+1) 
    fluxes in the y-direction. These are accessible via q.x(i,j) and q.y(i,j).
    
    \author Clancy Rowley
    \author $LastChangedBy$
    \date  4 Jul 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class Flux : public Field {
public:
    /// Constructor: allocate arrays based on the Grid dimensions
    Flux(const Grid& grid);

    /// Default constructor: set all dimensions to zero
    Flux();
    
    /// Constructor, making a copy of the data
    Flux(const Flux& q);

    /// Deallocate memory in the destructor
    ~Flux();

    /// Set all parameters and reallocate arrays based on the Grid dimensions
    void resize( const Grid& grid );
    
    /// Print the X and Y components to standard out (for debugging)
    void print();
    
    /// Copy assignment
    inline Flux& operator=(const Flux& q) {
        assert(q.Ngrid() == Ngrid());
        assert(q.Nx() == Nx());
        assert(q.Ny() == Ny());
        _data = q._data;
        return *this;
    }

    /// Copy assignment from double
    inline Flux& operator=(double a) {
        _data = a;
        return *this;
    }

    /// q(dir,i,j) refers to the flux in direction dir (X or Y) at edge (i,j)
    inline double& operator()(int lev, int dir, int i, int j) {
        assert( lev >= 0 && lev < Ngrid() );
        return _data(lev,getIndex(dir,i,j));
    };

    /// q(dir,i,j) refers to the flux in direction dir (X or Y) at edge (i,j)
    inline double operator()(int lev, int dir, int i, int j) const {
        assert( lev >= 0 && lev < Ngrid() );
        return _data(lev,getIndex(dir,i,j));
    };

    /// Type used for referencing elements
    typedef int index;
    
    /// f(ind) refers to the value corresponding to the given index ind
    inline double& operator()(int lev, index ind) {
        assert( lev >= 0 && lev < Ngrid() );
        assert( (ind >= 0) && (ind < _numFluxes) );
        return _data(lev, ind);
    }

    /// f(ind) refers to the value corresponding to the given index ind
    inline double operator()(int lev, index ind) const {
        assert( lev >= 0 && lev < Ngrid() );
        assert( (ind >= 0) && (ind < _numFluxes) );
        return _data(lev, ind);
    }

    /// q.x(dir,i) refers to the x-coordinate of the flux (dir,i,j)
    inline double x(int lev, int dir, int i) {
        assert( lev >= 0 && lev < Ngrid() );
        assert( (dir >= X) && (dir <= Y) );
        assert(i >= 0);
        if (dir == X) {
            assert( i < Nx()+1 );
            return getXEdge(lev,i);
        }
        else {
            assert( i < Nx() );
            return getXCenter(lev,i);
        }
    }

    /// q.x(ind) returns the x-coordinate of the flux specified by ind
    inline double x(int lev, index ind) {
        assert( lev >= 0 && lev < Ngrid() );
		int dir;
		if ( ind < _numXFluxes ) { 
			dir = 0; 
		} else {
			dir = 1;
		}
        int i = (ind - dir*_numXFluxes) / (Ny()+dir);
        return x(lev,dir,i);
    }
    
    /// q.y(ind) returns the x-coordinate of the flux specified by ind
    inline double y(int lev, index ind) {
        assert( lev >= 0 && lev < Ngrid() );
		int dir;
		if ( ind < _numXFluxes ) { 
			dir = 0; 
		} else {
			dir = 1;
		}
        int i = (ind - dir*_numXFluxes) / (Ny()+dir);
        int j = ind - dir*_numXFluxes - i * (Ny()+dir);
        return y(lev,dir,j);
    }
    
    /// q.y(dir,j) refers to the y-coordinate of the flux (dir,i,j)
    inline double y(int lev, int dir, int j) {
        assert( lev >= 0 && lev < Ngrid() );
        assert( (dir >= X) && (dir <= Y) );
        assert(j >= 0);
        if (dir == X) {
            assert( j < Ny() );
            return getYCenter(lev,j);
        }
        else {
            assert( j < Ny()+1 );
            return getYEdge(lev,j);
        }
    }

    /// Returns an index that refers to the first element
    inline index begin() const { return 0; }

    /// Returns an index that is one past the last element
    inline index end() const { return _numFluxes; }

    /// Returns an index for the first element in direction dir (X or Y)
    inline index begin(int dir) {
        assert ((dir >= X) && (dir <= Y));
        return dir * _numXFluxes;
    }

    /// Returns an index one past the last element in direction dir (X or Y)
    inline index end(int dir) {
        assert ((dir >= X) && (dir <= Y));
        if (dir == X) {
            return _numXFluxes;
        }
        else {
            return _numFluxes;
        }
    }
    
    /// Returns an index for the value in direction dir at point (i,j)
    inline index getIndex(int dir, int i, int j) const {
        assert(dir>=X  && dir<=Y);
        assert(i >= 0 && j >= 0);
        assert((dir == X) ? i < Nx()+1 : i < Nx());
        assert((dir == Y) ? j < Ny()+1 : j < Ny());
        // Tricky expression:
        //   j in [0..ny-1] for X fluxes (dir = X)
        //   j in [0..ny] for Y fluxes   (dir = Y)
        return dir * _numXFluxes + i * (Ny()+dir) + j;
    }
    
    /// f += g
    inline Flux& operator+=(const Flux& f) {
        assert(f.Ngrid() == Ngrid());
        assert(f.Nx() == Nx());
        assert(f.Ny() == Ny());
        for (unsigned int i=0; i < f._data.Size(); ++i) {
            _data(i) += f._data(i);
        }
        return *this;
    }

    /// f += a
    inline Flux& operator+=(double a) {
        for (unsigned int i=0; i < _data.Size(); ++i) {
            _data(i) += a;
        }
        return *this;
    }
    
    /// f -= g
    inline Flux& operator-=(const Flux& f) {
        assert(f.Ngrid() == Ngrid());
        assert(f.Nx() == Nx());
        assert(f.Ny() == Ny());
        for (unsigned int i=0; i < f._data.Size(); ++i) {
            _data(i) -= f._data(i);
        }
        return *this;
    }

    /// f -= a
    inline Flux& operator-=(double a) {
        for (unsigned int i=0; i < _data.Size(); ++i) {
            _data(i) -= a;
        }
        return *this;
    }
    
    /// f + g
    inline Flux operator+(const Flux& f) {
        assert(f.Ngrid() == Ngrid());
        assert(f.Nx() == Nx());
        assert(f.Ny() == Ny());
        Flux g = *this;
        g += f;
        return g;
    };
    
    /// f + a
    inline Flux operator+(double a) {
        Flux g = *this;
        g += a;
        return g;
    };

    /// f - g
    inline Flux operator-(const Flux& f) {
        assert(f.Ngrid() == Ngrid());
        assert(f.Nx() == Nx());
        assert(f.Ny() == Ny());
        Flux g = *this;
        g -= f;
        return g;
    };
    
    /// f - a
    inline Flux operator-(double a) {
        Flux g = *this;
        g -= a;
        return g;
    };

    /// f *= a
    inline Flux& operator*=(double a) {
        _data *= a;
        return *this;
    }

    /// f /= a
    inline Flux& operator/=(double a) {
        _data /= a;
        return *this;
    }

    /// f * a
    inline Flux operator*(double a) {
        Flux g = *this;
        g *= a;
        return g;
    }

    /// f / a
    inline Flux operator/(double a) {
        Flux g = *this;
        g /= a;
        return g;
    }

    /// Return Flux for a uniform flow with the specified magnitude and dir
    static Flux UniformFlow(
        const Grid& grid,
        double magnitude,
        double angle
    );

    /// Return Flux for moving flow given by element g\in TSE2 and a center of rotation (xCenter, yCenter)
    void setFlow(
        TangentSE2 g,
        double xCenter,
        double yCenter
    );

private:
    int _numXFluxes;
    int _numFluxes;
    Array::Array2<double> _data;
};  // class Flux

/// -f
inline Flux operator-(const Flux& f) {
    Flux g = f;
    g *= -1;
    return g;
}

/// a + f
inline Flux operator+(double a, const Flux& f) {
    Flux g = f;
    g += a;
    return g;
};
    
/// a - f
inline Flux operator-(double a, const Flux& f) {
    Flux g = -f;
    g += a;
    return g;
};

/// a * f
inline Flux operator*(double a, const Flux& f) {
    Flux g = f;
    g *= a;
    return g;
}

} // namespace ibpm

#endif /* _FLUX_H_ */
