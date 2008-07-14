#ifndef _FLUX_H_
#define _FLUX_H_

#include "Grid.h"
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
public:
    enum Dimension {X, Y};

	/// Allocate memory in the constructor
	Flux(const Grid& grid) :
		_grid(grid),
		_nx(grid.getNx()),
		_ny(grid.getNy()),
		_xdata(_nx+1,_ny),
		_ydata(_nx,_ny+1) {};
	
	/// Allocate new array, copy the data
	inline Flux(const Flux& q) :
		_grid(q._grid),
		_nx(q._nx),
		_ny(q._ny),
		_xdata(_nx+1,_ny),
		_ydata(_nx,_ny+1)
	{
		// copy data
		this->_xdata = q._xdata;
		this->_ydata = q._ydata;
	}

       	const Grid& getGrid() const{
		return _grid;
	}
	
	const Array<double, 2>& getDataX() const {
		return _xdata;
	}
	
	const Array<double, 2>& getDataY() const {
		return _ydata;
		}
	
    typedef Array<double,2>::iterator iterator;

    inline iterator begin(int dim) {
        assert(dim <= Y);
        switch (dim) {
            case X: return _xdata.begin();
            case Y: return _ydata.begin();
        }
    }
    
    inline iterator end(int dim) {
        assert(dim <= Y);
        switch (dim) {
            case X: return _xdata.end();
            case Y: return _ydata.end();
        }
    }
    
	/// Copy assignment
	inline Flux& operator=(const Flux& q) {
		assert(q._nx == this->_nx);
		assert(q._ny == this->_ny);
		this->_xdata = q._xdata;
		this->_ydata = q._ydata;
		return *this;
	}

	/// Copy assignment from double
	inline Flux& operator=(double a) {
		this->_xdata = a;
		this->_ydata = a;
		return *this;
	}

	/// q.x(i,j) refers to the x-flux at the left edge of cell (i,j)
	inline double& x(int i, int j) {
		assert(i < _nx+1);
		assert(j < _ny);
		return _xdata(i,j);
	}

	/// q.y(i,j) refers to the y-flux at the bottom edge of cell (i,j)
	inline double& y(int i, int j) {
		assert(i < _nx);
		assert(j < _ny+1);
		return _ydata(i,j);
	}

	/// Deallocate memory in the destructor
	~Flux() {};  // deallocation automatic for Blitz++ arrays?

	/// f += g
	inline Flux& operator+=(const Flux& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		this->_xdata += f._xdata;
		this->_ydata += f._ydata;
		return *this;
	}

	/// f += a
	inline Flux& operator+=(double a) {
		this->_xdata += a;
		this->_ydata += a;
		return *this;
	}
	
	/// f -= g
	inline Flux& operator-=(const Flux& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		this->_xdata -= f._xdata;
		this->_ydata -= f._ydata;
		return *this;
	}

	/// f -= a
	inline Flux& operator-=(double a) {
		this->_xdata -= a;
		this->_ydata -= a;
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
		this->_xdata *= a;
		this->_ydata *= a;
		return *this;
	}

	/// f /= a
	inline Flux& operator/=(double a) {
		this->_xdata /= a;
		this->_ydata /= a;
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

	/// Set *this to the curl of the argument, returning *this.
	Flux& curl(const Scalar& f);

	/// Set *this to the gradient of the argument, returning *this.
	Flux& gradient(const Scalar& f);
	
	/// Return the inner product of *this and the argument
	double dot(const Flux& q) const;
	
private:
	const Grid& _grid;
	const int _nx;
	const int _ny;
	Array<double,2> _xdata;
	Array<double,2> _ydata;	
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
