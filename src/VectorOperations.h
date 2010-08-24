#ifndef _VECTOROPERATIONS_H_
#define _VECTOROPERATIONS_H_

#include "Array.h"
using Array::Array2;

namespace ibpm {

class BC;
class Scalar;
class Flux;
class BoundaryVector;
class NavierStokesModel;

/*!
    \file VectorOperations.h
    \module VectorOperations

    \brief Collection of vector operations on Scalar and Flux objects.
    
    Operations include curl on scalars and fluxes, inner products of two 
    scalars or two fluxes, sine transform of scalars, and 'cross product'
    between fluxes/scalars.

    \author Clancy Rowley
    \author Zhanhua Ma
    \author $LastChangedBy$
    \date  15 Jul 2008
    \date $LastChangedDate$
    \version $Revision$
*/

/*! \brief Return the curl of Flux q, as a Scalar object.

The curl is defined only at the interior nodes, and this routine returns zero at the boundary nodes.
*/
Scalar Curl(const Flux& q);
void Curl(const Flux& q, Scalar& omega );
    
/// \brief Return the curl of Scalar f, as a Flux object. 
Flux Curl(const Scalar& f);
void Curl(const Scalar& f, Flux& q);

/// \brief Return the inner product of Scalar f and Scalar g, calculated over the finest grid only.
double FineGridInnerProduct( const Scalar& f, const Scalar& g );  
    
/// \brief Return the inner product of Scalar f and Scalar g.
double InnerProduct( const Scalar& f, const Scalar& g );  

/// \brief Return the inner product of Flux p and Flux q, calculated over the finest grid only.
double FineGridInnerProduct( const Flux& p, const Flux& q );

/// \brief Return the inner product of Flux p and Flux q.
double InnerProduct( const Flux& p, const Flux& q );
    
    
/// \brief Return the energy-equivalent inner product of two (Scalar) vorticity fields, calculated over the fine grid only.
/*  This is distinct from the standard L2 inner product so that the product of
 vorticity fields is equal to the product of fluxes (ie both yield KE)
 */
double FineGridVorticityInnerProduct( const Scalar& omega1, const Scalar& omega2, const NavierStokesModel& model );
    
/// \brief Return the energy-equivalent inner product of two (Scalar) vorticity fields.
double VorticityInnerProduct( const Scalar& omega1, const Scalar& omega2, const NavierStokesModel& model );

/// \brief Compute the Laplacian of f
void Laplacian( const Scalar& f, Scalar& g );
Scalar Laplacian( const Scalar& f );
void Laplacian(
    const Array2<double>& f,
    double dx,
    const BC& bc,
    Array2<double>& g );
        

/*! \brief Return the cross product of a Flux q and a Scalar f, as a Flux.

    q x f = (f v, -f u)
    
    This routine is designed so that the following identity holds at the
    discrete level:
        < a, q1 x q2 > = < q1, q2 x a >
*/
Flux CrossProduct(const Flux& q, const Scalar& f);

/*! \brief Return the cross product of two Flux objects, q1, q2, as a Scalar.

    q1 x q2 = u1 v2 - u2 v1

    This routine is designed so that the following identity holds at the
    discrete level:
        < a, q1 x q2 > = < q1, q2 x a >
*/
Scalar CrossProduct(const Flux& q, const Flux& p);

/// \brief Convert x-fluxes through edges to velocities at vertices
void FluxToXVelocity(const Flux& q, Scalar& u);

/// \brief Convert y-fluxes through edges to velocities at vertices
void FluxToYVelocity(const Flux& q, Scalar& v);

/// \brief Convert u-velocities at vertices to x-fluxes through edges.
/// Does not touch the y-component of the Flux q passed in.
void XVelocityToFlux(const Scalar& u, Flux& q);

/// \brief Convert v-velocities at vertices to y-fluxes through edges.
/// Does not touch the x-component of the Flux q passed in.
void YVelocityToFlux(const Scalar& v, Flux& q);

/// \brief Convert u- and v-velocities at vertices to fluxes through edges
void VelocityToFlux(const Scalar& u, const Scalar& v, Flux& q);

/// \brief Convert fluxes through edges to u- and v-velocities at vertices
void FluxToVelocity(const Flux& q, Scalar& u, Scalar& v);

} // namespace ibpm

#endif /* _VECTOROPERATIONS_H_ */
