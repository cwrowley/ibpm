#ifndef _FLUX_H_
#define _FLUX_H_

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

class Scalar;
class Grid;
class Flux;

/*!
	\file Flux.h
	\class Flux

	\brief Store a 2D array of fluxes, located at cell edges

	\author Clancy Rowley
	\author $LastChangedBy: $
	\date  4 Jul 2008
	\date $LastChangedDate: $
	\version $Revision: $
*/

class Flux {
public:
	/// Allocate memory in the constructor
	Flux(const Grid& grid);

	/// Deallocate memory in the destructor
	~Flux();
	
	/// Return number of cells in x-direction
	inline int getNx() const { return _nx; }
	
	/// Return number of cells in y-direction
	inline int getNy() const { return _ny; }
	
	/// Set *this to the curl of the argument, returning *this.
	Flux curl(const Scalar& f);

	/// Set *this to the gradient of the argument, returning *this.
	Flux gradient(const Scalar& f);
	
	/// Return the inner product of *this and the argument
	double dot(const Flux& q) const;
	
	// TODO: overloaded operators for addition, etc
	
private:
	const Grid* _grid;
	int _nx;
	int _ny;
	Array<double,2> _x;
	Array<double,2> _y;
};

#endif /* _FLUX_H_ */
