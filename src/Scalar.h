#ifndef _SCALAR_H_
#define _SCALAR_H_

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

class Grid;
class Flux;

/*!
	\file Scalar.h
	\class Scalar

	\brief Store a 2D array of scalar values, located at cell centers.
	
	Also provides scalar-valued math operations, such as curl and divergence
	of a Flux, discrete sin transform, and Laplacian and its inverse, as well
	as inner product of scalar fields.

	\author Clancy Rowley
	\author $LastChangedBy: $
	\date  4 Jul 2008
	\date $LastChangedDate: $
	\version $Revision: $
*/

class Scalar {
public:
	/// Allocate memory for the 2D array
	Scalar(const Grid& grid);

	// Scalar f, f1, f2;
	// double a;
	
	// TO DO: overload operator f(i,j)
	
	// TODO: overload operator for f1 + f2
	
	// TODO: overload operator for f1 - f2

	// TODO: overload operator for a + f
	
	// TODO: overload operator for a - f

	// TODO: overload operator for f + a
	
	// TODO: overload operator for f - a

	// TODO: overload operator for a * f
	
	// TODO: overload operator for f * a

	// TODO: overload operator for f / a
	
	/// Set *this to the curl of q
	Scalar curl(const Flux& q);
	
	/// set *this to the discrete sin transform of f
	Scalar sinTransform(const Scalar& f);

	/// set *this to the inverse discrete sin transform of f
	Scalar sinTransformInv(const Scalar& f);
	
	/// set *this to the Laplacian of f
	Scalar laplacian(const Scalar& f);

	/// set *this to the inverse Laplacian of f
	Scalar laplacianInverse(const Scalar& f);

	/// set *this to the divergence of q
	Scalar divergence(const Flux& q);

	/// return the inner product of f and *this
	double dot(const Scalar& f);
	
private:
	
	const Grid* _grid;
	int _nx;
	int _ny;
	Array<double,2> _data;
	
	// Declare variables for implementation here
};

#endif /* _SCALAR_H_ */
