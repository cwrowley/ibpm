#ifndef _VECTOROPERATIONS_H_
#define _VECTOROPERATIONS_H_

class Scalar;
class Flux;
class BoundaryVector;

/*!
	\file VectorOperations.h
	\module VectorOperations

	\brief Do vector operations on Scalar and Flux objects.
	
	Operations include curl on scalars and fluxes, inner products of two scalars or two fluxes,
	sine transform of scalars, and 'cross product' between fluxes/scalars.  
		
	Functions in this module are friends of classes Scalar and Flux.	 

	\author Clancy Rowley
	\author $LastChangedBy: zma $
	\date  15 Jul 2008
	\date $LastChangedDate: $
	\version $Revision: $
*/

/*! \brief Return the curl of Flux q, as a Scalar object.

The curl is defined only at the interior nodes, and this routine returns zero at the boundary nodes.
*/
Scalar curl(const Flux& q); 		 				 
	
/// Return the curl of Scalar f, as a Flux object. 
Flux curl(const Scalar& f);
	
/// Return the inner product of Scalar f and Scalar g.
double InnerProduct( const Scalar& f, const Scalar& g );

/// Return the inner product of Flux p and Flux q.
double InnerProduct( const Flux& p, const Flux& q );

/*! Return the sine transform of a Scalar object using DST-I.
(fftw library is used (real fft kind:RODFT00); only interior nodes are considered.)
*/
Scalar sinTransform(const Scalar& f);

/// Return the inverse sine transform of a Scalar object using fft.(fftw library
/// is used; only interior nodes are considered.)
//Scalar sinTransformInv(const Scalar& f);

/// Return the inner product of BoundaryVectors x and y.
// double InnerProduct( const BoundaryVector& x, const BoundaryVector& y );

/// Return the average of a flux in x direction, taken at nodes, as a Scalar object. 
/// BCs: By supposing zero 'ghost flux.x' outside of the domain.  
Scalar fluxXAverage(const Flux& q);

/// Return the average of a flux in y direction, taken at nodes, as a Scalar object. 
/// BCs: By supposing zero 'ghost flux.y' outside of the domain.  
Scalar fluxYAverage(const Flux& q);

/*! \brief Return the 'cross product'  of a Flux object q and a Scalar object f, as a Flux object.

	Step 1: Take the averages of q.x and q.y, respectively, at nodes.	
	Step 2: Cross product of (q.xavg, q.yavg, 0) and (0, 0, f).
	Step 3: Take the averages of the cross product at edges and return the flux. 
	
	BCs: The returned flux.X is zero at left and right boundaries (i = 0, nx);
	the returned flux.Y is zero at top and bottom boundaries (j = 0, ny).
*/
Flux crossproduct(const Flux& q, const Scalar& f);

/*! \brief Return the 'cross product' q cross p of two Flux objects, q, p, as a Scalar object.

	Step 1: Take the averages of q.x, q.y, p.x, p.y, respectively, at nodes.
	Step 2: Cross product of (q.xavg, q.yavg, 0) and (p.xavg, p.yavg, 0), returns a Scalar.
	
	BCs: Entries at the four boundary edges are computed with the assumption
	that the related 'ghost fluxes' outside of the domain are zero.   
*/
Scalar crossproduct(const Flux& q, const Flux& p);

#endif /* _VECTOROPERATIONS_H_ */
