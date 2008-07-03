/*! \file NavierStokesModel.h */

#ifndef _NAVIERSTOKES_H_
#define _NAVIERSTOKES_H_

class Geometry;

/*! \brief Define operators for the Navier-Stokes equations.
 *  
 */
class NavierStokesModel {
private:
	Geometry* _geometry;
};

//! Full nonlinear Navier-Stokes equations.
class NonlinearNavierStokes : public NavierStokesModel {};

//! Navier-Stokes equations linearized about an equilibrium point.
class LinearizedNavierStokes : public NavierStokesModel {};

//! Adjoint Navier-Stokes equations, linearized about an equilibrium point.
class AdjointNavierStokes : public NavierStokesModel {};


#endif /* _NAVIERSTOKES_H_ */
