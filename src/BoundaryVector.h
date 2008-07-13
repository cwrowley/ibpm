#ifndef _BOUNDARYVECTOR_H_
#define _BOUNDARYVECTOR_H_

/*!
	\file BoundaryVector.h
	\class BoundaryVector

	\brief Store vectors located at boundary points
	
	Examples are forces, velocities, and coordinates of points on a body.

	\author Clancy Rowley
	\author $LastChangedBy: $
	\date  7 Jul 2008
	\date $LastChangedDate: $
	\version $Revision: $
*/


class BoundaryVector {
public:
    enum Dimension {X, Y};

	/// \brief Constructor, allocating memory for a body with
	/// n points on the boundary.
	BoundaryVector(int n) :
	    _numPoints(n),
	    _data(_numPoints, 2),
	    _xdata( _data(Range::all(), 0) ),
	    _ydata( _data(Range::all(), 1) ) {};

	/*! \brief Constructor using pre-existing data, as a 1d array.

	Here, n is the number of boundary points, and *data has size 2n.
	*/
	BoundaryVector(int n, double* data);
	// TODO: Implement this, and write tests

	/// Return the number of boundary points
	inline int getNumPoints() { return _numPoints; }
	
	/*! \brief Return the number of elements in the array.
	
	(twice the number of boundary points)
	*/
	inline int getSize() { return 2*_numPoints; }
	
	/// Reference to the x-component at boundary point i
	inline double& x(int i) { return _data(i,0); }
	
	/// Reference to the y-component at boundary point i
	inline double& y(int i) { return _data(i,1); }
	
	/// Return a pointer to the data, expressed as a C-style array.
	double* flatten();
	// TODO: Implement this, and write tests
	
    typedef Array<double,1>::iterator iterator;

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
    
	/// Return the dot product of *this and the argument
	double dot(BoundaryVector& f);
	// TODO: Implement this, and write tests

	/// Copy assignment
	inline BoundaryVector& operator=(const BoundaryVector& f) {
		assert(f._numPoints == this->_numPoints);
		this->_data = f._data;
		return *this;
	};

	/// Copy assignment from double
	inline BoundaryVector& operator=(double a) {
		this->_data = a;
		return *this;
	};

	/// f += g
	inline BoundaryVector& operator+=(const BoundaryVector& f) {
		assert(f._numPoints == this->_numPoints);
		this->_data += f._data;
		return *this;
	}

	/// f -= g
	inline BoundaryVector& operator-=(const BoundaryVector& f) {
		assert(f._numPoints == this->_numPoints);
		this->_data -= f._data;
		return *this;
	}

	/// f *= a
	inline BoundaryVector& operator*=(double a) {
		this->_data *= a;
		return *this;
	}
	
	/// f /= a
	inline BoundaryVector& operator/=(double a) {
		this->_data /= a;
		return *this;
	}
	
	/// f + g
	inline BoundaryVector operator+(const BoundaryVector& g) {
		BoundaryVector h = *this;
		h += g;
		return h;
	}
	
	/// f - g
	inline BoundaryVector operator-(const BoundaryVector& g) {
		BoundaryVector h = *this;
		h -= g;
		return h;
	}
	
	/// f * a
	inline BoundaryVector operator*(double a) {
		BoundaryVector g = *this;
		g *= a;
		return g;
	}
	
	/// f / a
	inline BoundaryVector operator/(double a) {
		BoundaryVector g = *this;
		g /= a;
		return g;
	}

private:
	int _numPoints;
    Array<double,2> _data;
    Array<double,1> _xdata; // 1d slice of all x components of _data
    Array<double,1> _ydata; // 1d slice of all y components of _data
};

/// -f
inline BoundaryVector operator-(const BoundaryVector& f) {
	BoundaryVector g = f;
	g *= -1;
	return g;
}

/// a * f
inline BoundaryVector operator*(double a, const BoundaryVector& f) {
	BoundaryVector g = f;
	g *= a;
	return g;
}
	



#endif /* _BOUNDARYVECTOR_H_ */
