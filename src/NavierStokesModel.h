#ifndef _NAVIERSTOKES_H_
#define _NAVIERSTOKES_H_

#include "Grid.h"
#include "Geometry.h"
#include "Scalar.h"
#include "Flux.h"
#include "BaseFlow.h"
#include "BoundaryVector.h"
#include "State.h"
#include "VectorOperations.h"
#include "Regularizer.h"
#include "EllipticSolver.h"
#include <math.h>

namespace ibpm {

/*!
\file NavierStokesModel.h
\class NavierStokesModel

\brief Define operators for the Navier-Stokes equations

Note: This is an abstract base class, and may not be instantiated directly

\author Clancy Rowley
\date  3 Jul 2008

$Revision$
$LastChangedDate$
$LastChangedBy$
$HeadURL$
*/

class NavierStokesModel {
public:
    /*! \brief Constructor, given a specified grid, geometry,
        and potential flow flux
    */
    NavierStokesModel(
        const Grid& grid,
        const Geometry& geometry,
		double Reynolds,
        const BaseFlow& q_potential
    );

    /// \brief Constructor, assuming potential flow part is zero
    NavierStokesModel(
        const Grid& grid,
        const Geometry& geometry,
		double Reynolds
    );

    ~NavierStokesModel();

    /// Perform initial calculations needed to use model
    void init();

    /// \brief Return true if the geometry has moving bodies
    bool isTimeDependent() const;
    bool geTimeDependent() const;  // geometry TD?
    bool bfTimeDependent() const;  // baseflow TD?


    /// \brief Return the number of points in the geometry
    int getNumPoints() const;

    /// \brief Return the right-hand side b of the constraint equations.
    /// Here, this is the velocity of the bodies minus the base flow velocity
    BoundaryVector getConstraints() const;

    /// \brief Return the boundary points in the geometry
    BoundaryVector getPoints() const;

    /// \brief Update operators, for time-dependent models
    void updateOperators( double time );

//    /// Return a pointer to the associated Geometry
//    inline const Geometry& getGeometry() const { return _geometry; }
//
    /// Return a pointer to the associated Grid
    inline const Grid& getGrid() const { return _grid; }

    /// Compute omega = B(f) as in (14)
    void B(const BoundaryVector& f, Scalar& omega ) const;

    /// Compute f = C(omega) as in (14)
    void C(const Scalar& omega, BoundaryVector& f) const;

	/// \brief Return the constant alpha = 1/ReynoldsNumber
    double getAlpha() const;

    /// Return the angle of attack of the baseflow
    inline double getAlphaBF() const { return _baseFlow.getAlpha(); }

    /// Compute flux q from vorticity omega, including base flow q0
    void computeFlux(const Scalar& omega, Flux& q ) const;

    /// \brief Compute flux q from the vorticity omega, including base flow
    void refreshState( State& x ) const;

	/*! \brief Given the vorticity omega, return the streamfunction psi.

	 Assumes psi = 0 on the boundary, and does not add in potential flow solution
	 */
    Scalar vorticityToStreamfunction(const Scalar& omega) const;

    /// Return the BaseFlow
    inline const BaseFlow& getBaseFlow() const { return _baseFlow; }

    // Return flux corresponding to forces
    Flux force2flux(BoundaryVector& f);

    // Return boundary velocity corresponding to flux
    BoundaryVector flux2boundary(Flux &q);

private:
    // Perhaps copy constructor should be made private to prevent future memory deallocation errors

    BoundaryVector getBaseFlowBoundaryVelocities() const;
    void computeFluxWithoutBaseFlow(const Scalar& omega, Flux& q ) const;

    // data
    const Grid& _grid;
    const Geometry& _geometry;
    Regularizer _regularizer;
    BaseFlow _baseFlow;
	double _ReynoldsNumber;
    PoissonSolver _poisson;
    bool _hasBeenInitialized;
};

} // namespace ibpm

#endif /* _NAVIERSTOKES_H_ */
