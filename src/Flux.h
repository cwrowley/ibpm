#ifndef _FLUX_H_
#define _FLUX_H_

class Scalar;
class Grid;

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
};

#endif /* _FLUX_H_ */
