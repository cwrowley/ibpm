#ifndef _VECTOROPERATIONS_H_
#define _VECTOROPERATIONS_H_

class Scalar;
class Flux;

/*!
	\file VectorOperations.h
	\module VectorOperations

	\brief Do vector operations on Scalar and Flux objects.
	
	Operations include curl and divergence on scalars, and curl and gradient on
	fluxes.
		
	Functions in this module are friends of classes Scalar and Flux.
	 
	Blitz++ library is used for array operations.

	\author Clancy Rowley
	\author $LastChangedBy: zma $
	\date  15 Jul 2008
	\date $LastChangedDate: $
	\version $Revision: $
*/

/// Return the curl of Flux q, as a Scalar object. Return zero at boundary nodes.
Scalar curl(const Flux& q); 		 				 
	
/// Return the divergence of Flux q, as a Scalar object
Scalar divergence(const Flux& q);

/// Return the curl of Scalar f, as a Flux object. 
Flux curl(const Scalar& f);

/// Return the gradient of Scalar f, as a Scalar object.
Flux gradient(const Scalar& f);


#endif /* _VECTOROPERATIONS_H_ */
