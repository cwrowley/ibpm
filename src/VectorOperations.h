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

/*! \brief Return the curl of Flux q, as a Scalar object.

The curl is defined only at the interior nodes, and this routine returns zero at the boundary nodes.
*/
Scalar curl(const Flux& q); 		 				 
	
/// Return the curl of Scalar f, as a Flux object. 
Flux curl(const Scalar& f);

#endif /* _VECTOROPERATIONS_H_ */
