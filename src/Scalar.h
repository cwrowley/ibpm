#ifndef _SCALAR_H_
#define _SCALAR_H_

#include <blitz/array.h>
#include "Grid.h"
#include "Flux.h"

BZ_USING_NAMESPACE(blitz)

//class Flux;

/*!
	\file Scalar.h
	\class Scalar

	\brief Store a 2D array of scalar values, located at cell nodes.
	
	For a grid with nx cells in the x-direction and ny cells in the y-direction,
	there are (nx+1)*(ny+1) nodes.  
	
	There are (nx-1)*(ny-1) inner nodes, and 2*(nx+ny) boundary nodes.
	
	Also provides scalar-valued math operations, such as discrete sin transform, 
	and Laplacian and its inverse, as well as inner product of scalar fields.

	\author Clancy Rowley
	\author $LastChangedBy$
	\date  4 Jul 2008
	\date $LastChangedDate$
	\version $Revision$
*/

class Scalar {
	/// give VectorOperation functions access to private data
    friend Scalar curl(const Flux& q);
    friend Scalar divergence(const Flux& q);
	friend Flux curl(const Scalar& f);
	friend Flux gradient(const Scalar& f); 

public:
	/// Allocate memory for the 2D array
	Scalar(const Grid& grid) :
		_grid(grid),
		_nx(grid.getNx()),
		_ny(grid.getNy()),
		_data(_nx+1,_ny+1) {};
	
	/// Allocate new array, copy the data
	inline Scalar(const Scalar& f) :
		_grid(f._grid),
		_nx(f._nx),
		_ny(f._ny),
		_data(_nx+1,_ny+1)
	{
		// copy data
		this->_data = f._data;
	}
	
	/// Deallocate memory in the destructor
	~Scalar() {
		// deallocation automatic for Blitz++ arrays?
	};

	const Grid& getGrid() const{
		return _grid;
	}
	
    // const Array<double, 2>& getData() const {
    //  return _data;
    // }
	
    typedef Array<double,2>::iterator iterator;

    inline iterator begin() {
        return _data.begin();
    }
    
    inline iterator end() {
        return _data.end();
    }
    
	/// Copy assignment
	inline Scalar& operator=(const Scalar& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		this->_data = f._data;
		return *this;
	};

	/// Copy assignment from double
	inline Scalar& operator=(double a) {
		this->_data = a;
		return *this;
	};

	/// f(i,j) refers to the value at index (i,j)
	inline double& operator()(int i, int j) {
		assert(i>=0  && i<=_nx);
		assert(j>=0  && j<=_ny);
		return _data(i,j);
	};
	
	/// f += g
	inline Scalar& operator+=(const Scalar& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		this->_data += f._data;
		return *this;
	}

	/// f += a
	inline Scalar& operator+=(double a) {
		this->_data += a;
		return *this;
	}

	/// f -= g
	inline Scalar& operator-=(const Scalar& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		this->_data -= f._data;
		return *this;
	}

	/// f -= a
	inline Scalar& operator-=(double a) {
		this->_data -= a;
		return *this;
	}

	/// f + g
	inline Scalar operator+(const Scalar& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		Scalar g = *this;
		g += f;
		return g;
	};
	
	/// f + a
	inline Scalar operator+(double a) {
		Scalar g = *this;
		g += a;
		return g;
	};
	
	/// f - g
	inline Scalar operator-(Scalar& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		Scalar g = *this;
		g -= f;
		return g;
	};
	
	/// f - a
	inline Scalar operator-(double a) {
		Scalar g = *this;
		g -= a;
		return g;
	};
	
	/// f *= g
	inline Scalar& operator*=(const Scalar& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		this->_data *= f._data;
		return *this;
	}

	/// f *= a
	inline Scalar& operator*=(double a) {
		this->_data *= a;
		return *this;
	}

	/// f * g
	inline Scalar operator*(const Scalar& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		Scalar g = *this;
		g *= f;
		return g;
	};
	
	/// f * a
	inline Scalar operator*(double a) {
		Scalar g = *this;
		g *= a;
		return g;
	};
	
	/// f /= g
	inline Scalar& operator/=(const Scalar& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		this->_data /= f._data;
		return *this;
	}

	/// f /= a
	inline Scalar& operator/=(double a) {
		this->_data /= a;
		return *this;
	}

	/// f / g
	inline Scalar operator/(Scalar& f) {
		assert(f._nx == this->_nx);
		assert(f._ny == this->_ny);
		Scalar g = *this;
		g /= f;
		return g;
	};
	
	/// f / a
	inline Scalar operator/(double a) {
		Scalar g = *this;
		g /= a;
		return g;
	};
	
	
	/// set *this to the discrete sin transform of f
	Scalar& sinTransform(const Scalar& f);

	/// set *this to the inverse discrete sin transform of f
	Scalar& sinTransformInv(const Scalar& f);
	
	/// set *this to the Laplacian of f
	Scalar& laplacian(const Scalar& f);

	/// set *this to the inverse Laplacian of f
	Scalar& laplacianInverse(const Scalar& f);

	
	/// return the inner product of f and *this
	double dot(const Scalar& f);
    
private:
	
	const Grid& _grid;
	const int _nx;
	const int _ny;
	Array<double,2> _data;
	
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
};
	
/// a - f
inline Scalar operator-(double a, Scalar& f) {
	Scalar g = -f;
	g += a;
	return g;
};
	
/// a * f
inline Scalar operator*(double a, Scalar& f) {
	Scalar g = f;
	g *= a;
	return g;
};

/// a / f
inline Scalar operator/(double a, Scalar& f) {
	Scalar g(f.getGrid());
	g = a;
	g /= f;
	return g;
};

#endif /* _SCALAR_H_ */
