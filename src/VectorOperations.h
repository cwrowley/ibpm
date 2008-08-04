#ifndef _VECTOROPERATIONS_H_
#define _VECTOROPERATIONS_H_

class Scalar;
class Flux;

/*!
	\file VectorOperations.h
	\module VectorOperations

	\brief Do vector operations on Scalar and Flux objects.
	
	Operations include curl and divergence on scalars, curl and gradient on
	fluxes, as well as inner products of two scalars or two fluxes.  
		
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
double InnerProduct (const Scalar& f, const Scalar& g);

/// Return the inner product of Flux p and Flux q.
double InnerProduct (const Flux& p, const Flux& q);

#endif /* _VECTOROPERATIONS_H_ */
