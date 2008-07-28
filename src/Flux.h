#ifndef _FLUX_H_
#define _FLUX_H_

#include "Grid.h"
#include "Direction.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

class Scalar;

/*!
\file Flux.h
\class Flux

\brief Store a 2D array of fluxes, located at cell edges

For a grid with nx cells in the x-direction and ny cells in the y-direction,
there are (nx+1,ny) fluxes in the x-direction, and (nx,ny+1) fluxes in the
y-direction. These are accessible via q.x(i,j) and q.y(i,j).

\author Clancy Rowley
\author $LastChangedBy: $
\date  4 Jul 2008
\date $LastChangedDate: $
\version $Revision: $
*/

class Flux {
	/// give VectorOperation functions access to private data
    friend Scalar curl(const Flux& q);
    friend Scalar divergence(const Flux& q);
	friend Flux curl(const Scalar& f);
	friend Flux gradient(const Scalar& f); 
	
public:
	/// Allocate memory in the constructor
	Flux(const Grid& grid) :
		_grid(grid),
		_nx(grid.getNx()),
		_ny(grid.getNy()),
		_numXFluxes(_nx * _ny + _ny),
		_numFluxes(2 * _nx * _ny + _nx + _ny),
        _data(_numFluxes) {};
	
	/// Allocate new array, copy the data
	inline Flux(const Flux& q) :
		_grid(q._grid),
		_nx(q._nx),
		_ny(q._ny),
		_numXFluxes(q._numXFluxes),
		_numFluxes(q._numFluxes),
		_data(_numFluxes)
	{
		// copy data
		this->_data = q._data;
	}

       	const Grid& getGrid() const{
		return _grid;
	}
	
	/// Deallocate memory in the destructor
	~Flux() {};  // deallocation automatic for Blitz++ arrays?

    // typedef Array<double,2>::iterator iterator;
    // 
    // inline iterator begin(int dim) {
    //     assert(dim >= X && dim <= Y);
    //     switch (dim) {
    //         case X: return _xdata.begin();
    //         case Y: return _ydata.begin();
    //     }
    // }
    // 
    // inline iterator end(int dim) {
    //     assert(dim <= Y);
    //     switch (dim) {
    //         case X: return _xdata.end();
    //         case Y: return _ydata.end();
    //     }
    // }
    
	/// Copy assignment
	inline Flux& operator=(const Flux& q) {
		assert(q._nx == this->_nx);
		assert(q._ny == this->_ny);
		this->_data = q._data;
		return *this;
	}

	/// Copy assignment from double
	inline Flux& operator=(double a) {
		this->_data = a;
		return *this;
	}

    /// q.x(dir,i) refers to the x-coordinate of the flux (dir,i,j)
    inline double x(int dir, int i) {
        assert( (dir >= X) && (dir <= Y) );
        assert(i >= 0);
        if (dir == X) {
            assert( i < _nx+1 );
            return _grid.getXEdge(i);
        } else {
            assert( i < _nx );
            return _grid.getXCenter(i);
        }
    }

    /// q.x(ind) returns the x-coordinate of the flux specified by ind
    // inline double x(index ind) {
    //     
    // }
    
    /// q.y(dir,j) refers to the y-coordinate of the flux (dir,i,j)
    inline double y(int dir, int j) {
        assert( (dir >= X) && (dir <= Y) );
        assert(j >= 0);
        if (dir == X) {
            assert( j < _ny );
            return _grid.getYCenter(j);
        } else {
            assert( j < _ny+1 );
            return _grid.getYEdge(j);
        }
    }

	/// Type used for referencing elements
    typedef int index;
    
    /// f(ind) refers to the value corresponding to the given index ind
	inline double& operator()(index ind) {
        assert( (ind >= 0) && (ind < _numFluxes) );
        return _data(ind);
	}

    /// Returns an index that refers to the first element
    inline index begin() { return 0; }

    /// Returns an index that is one past the last element
    inline index end() { return _numFluxes; }

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
        } else {
            return _numFluxes;
        }
    }
    
    /// Returns an index for the value in direction dir at point (i,j)
	inline index getIndex(int dir, int i, int j) const {
	    assert(dir>=X  && dir<=Y);
        assert(i >= 0 && j >= 0);
        assert((dir == X) ? i < _nx+1 : i < _nx);
        assert((dir == Y) ? j < _ny+1 : j < _ny);
        // Tricky expression:
        //   j in [0..ny-1] for X fluxes (dir = X)
        //   j in [0..ny] for Y fluxes   (dir = Y)
		return dir * _numXFluxes + i * (_ny+dir) + j;
	}
	
    /// q(dir,i,j) refers to the flux in direction dir (X or Y) at edge (i,j)
	inline double& operator()(int dir, int i, int j) {
        return _data(this->getIndex(dir,i,j));
	};

    /// q(dir,i,j) refers to the flux in direction dir (X or Y) at edge (i,j)
	inline double operator()(int dir, int i, int j) const {
        return _data(this->getIndex(dir,i,j));
	};

	/// f += g
	inline Flux& operator+=(const Flux& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		this->_data += f._data;
		return *this;
	}

	/// f += a
	inline Flux& operator+=(double a) {
		this->_data += a;
		return *this;
	}
	
	/// f -= g
	inline Flux& operator-=(const Flux& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		this->_data -= f._data;
		return *this;
	}

	/// f -= a
	inline Flux& operator-=(double a) {
		this->_data -= a;
		return *this;
	}
	
	/// f + g
	inline Flux operator+(const Flux& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
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
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
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
		this->_data *= a;
		return *this;
	}

	/// f /= a
	inline Flux& operator/=(double a) {
		this->_data /= a;
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
	
	/// Return the inner product of *this and the argument
	double dot(const Flux& q) const;

	
private:
	const Grid& _grid;
	const int _nx;
	const int _ny;
    const int _numXFluxes;
    const int _numFluxes;
	Array<double,1> _data;
};

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


#endif /* _FLUX_H_ */
