#ifndef _REGULARIZER_H_
#define _REGULARIZER_H_

/*!
\file Regularizer.h
\class Regularizer

\brief Define regularization operations between grid and boundary data

Regularization (smearing) lifts values defined on the boundary to values
defined on a grid, and interpolation takes data defined on a grid and
interpolates it to the boundary.

Specifically, if \f$ u(x) \f$ are values of u defined on a grid, and if
\f$ u(\xi) \f$ denotes values on a boundary, then regularization defines

\f[
u(x) = \int_\Omega u(\xi) \delta(x-\xi)\,d\xi
\f]

where a discrete approximation of the delta function is used, and
interpolation defines

\f[
u(\xi) = \int_\Omega u(x) \delta(x-\xi)\,dx
\f]

These integrals are discretized using a regularized version of the delta
function, with finite support, as in (14) of Taira & Colonius (J Comput Phys,
2007).

\author Clancy Rowley
\author $LastChangedBy: $
\date  7 Jul 2008
\date $LastChangedDate: $
\version $Revision: $
*/

class Regularizer {
public:
	/// Constructor: allocates memory for internal data structures
	Regularizer(const Grid& grid, const Geometry& geometry);
	
	/// Destructor: deallocate memory for internal data structures
    ~Regularizer();
	
	/// Update operators, for instance when the position of the bodies changes
	void update();
	
	/// Smear boundary data to grid
	Flux toGrid(const BoundaryVector& u);

	/// Interpolate grid data to boundary
	BoundaryVector toBoundary(const Flux& u);

private:
	double deltaFunction(double r);
	const Grid& _grid;
	const Geometry& _geometry;
};

#endif /* _REGULARIZER_H_ */
