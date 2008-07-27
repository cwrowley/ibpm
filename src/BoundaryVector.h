#ifndef _BOUNDARYVECTOR_H_
#define _BOUNDARYVECTOR_H_

#include <blitz/array.h>
#include "Direction.h"

BZ_USING_NAMESPACE(blitz)

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

class BoundaryVectorIndex;

class BoundaryVector {
public:
	/// \brief Constructor, allocating memory for a body with
	/// n points on the boundary.
	BoundaryVector(int n) :
	    _numPoints(n),
	    _data( XY, _numPoints) {};

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
	
	/// f(dir,i) refers to the value in direction dir (X or Y) at point i
	inline double& operator()(int dir, int i) {
		assert(dir>=X  && dir<=Y);
		assert(i>=0  && i<_numPoints);
		return _data(dir,i);
	};
	
    typedef BoundaryVectorIndex index;

	inline double& operator()(index ind) {
        return _data(ind.getIndex());
	}

    inline index getIndex(int dir, int i) {
		assert(dir>=X  && dir<=Y);
		assert(i>=0  && i<_numPoints);
        index ind = getIndex(*this, dir, i);
        return ind;
    }
	/// Return a pointer to the data, expressed as a C-style array.
	double* flatten();
	// TODO: Implement this, and write tests
	
    // typedef Array<double,1>::iterator iterator;
    // 
    // inline iterator begin(int dim) {
    //     assert(dim <= Y);
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
	
// /// Type for referencing an element of a BoundaryVector
// class BoundaryVectorIndex {
// public:
//     BoundaryVectorIndex() {
//         _numPoints = 0;
//         ind = 0,0;
//     }
// 
//     BoundaryVectorIndex(BoundaryVector& f) {
//         _numPoints = f.getNumPoints();
//         ind = 0,0;
//     }
// 
//     // public data: should clean this up later
//     TinyVector<2,int> ind;
// 
//     inline BoundaryVectorIndex operator++(BoundaryVectorIndex& index) {
//         if (++(index.ind(1)) >= _numPoints) {
//             index.ind(1) = 0;
//             ++(index.ind(0));
//         };
//     }
//     
// private:
//     int _numPoints
// }
// 
// /// Return an index for a given pair (dir, i), for a specified BV f
// BoundaryVectorIndex getIndex(BoundaryVector& f, Direction dir, int i) {
//     index = BoundaryVectorIndex(f);
//     index.ind = dir, i;
//     return index;
// }


#endif /* _BOUNDARYVECTOR_H_ */
