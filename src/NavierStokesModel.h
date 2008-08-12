#ifndef _NAVIERSTOKES_H_
#define _NAVIERSTOKES_H_

#include "Grid.h"

class Scalar;
class Geometry;
class BoundaryVector;
class State;
class Flux;

/*!
\file NavierStokesModel.h
\class NavierStokesModel

\brief Define operators for the Navier-Stokes equations

\author Clancy Rowley
\date  3 Jul 2008

$Revision$
$LastChangedDate$
$LastChangedBy$
$HeadURL$
*/

class NavierStokesModel {
public:
	/*! \brief Construcor, given a specied grid, geometry, and Reynolds number
	*/
	NavierStokesModel(
		const Grid& grid,
		const Geometry& geometry,
		double Reynolds
	);

	~NavierStokesModel() {}

	/// Return a pointer to the associated Geometry
	const Geometry* getGeometry() const { return &_geometry; }

	/// Return a pointer to the associated Grid
	const Grid* getGrid() const { return &_grid; }

	/// Return a pointer to the eigenvalues of the linear term L
	Scalar* getLambda() const { return _lambda; }
	
	/// Transform to eigenvectors of L (discrete sin transform)
	Scalar S(const Scalar& g) const;

	/*! \brief Inverse transform of S.
	Note that for the discrete sin transform, the S^{-1} should equal S, but
	in some implementations (e.g. FFTW), these may differ by a normalization
	constant.
	*/
	Scalar Sinv(const Scalar& ghat) const;
	
	/// Compute gamma = B(f) as in (14)
	Scalar B(const BoundaryVector& f) const;
	
	/// Compute f = C(gamma) as in (14)
	BoundaryVector C(const Scalar& gamma) const;
	
	/*! \brief Compute nonlinear terms y = N(x)
	Pure virtual function: must be overridden by subclasses.
	*/
	virtual Scalar nonlinear(const State& x) const = 0;
	
	/// Compute flux q from circulation gamma
    void computeFlux( const Scalar& gamma, Flux& q ) const;

protected:
	/*! \brief Compute bilinear term, used by subclasses
	(e.g. y = N(q) = bilinear(q,q) )
	*/
	Scalar bilinear(const State& x1, const State& x2);
	
private:
	const Geometry& _geometry;
	const Grid& _grid;
	Scalar* _lambda;
};

//! Full nonlinear Navier-Stokes equations.
class NonlinearNavierStokes : public NavierStokesModel {};

//! Navier-Stokes equations linearized about an equilibrium point.
class LinearizedNavierStokes : public NavierStokesModel {};

//! Adjoint Navier-Stokes equations, linearized about an equilibrium point.
class AdjointNavierStokes : public NavierStokesModel {};


#endif /* _NAVIERSTOKES_H_ */
