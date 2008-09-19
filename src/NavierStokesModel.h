#ifndef _NAVIERSTOKES_H_
#define _NAVIERSTOKES_H_

#include "Grid.h"
#include "Geometry.h"
#include "Scalar.h"
#include "Flux.h"
#include "BoundaryVector.h"
#include "State.h"
#include "VectorOperations.h"
#include "Regularizer.h"
#include "EllipticSolver.h"
#include "Model.h"
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

class NavierStokesModel : public Model {
public:
    /*! \brief Constructor, given a specified grid, geometry, Reynolds number,
        and potential flow flux
    */
    NavierStokesModel(
        const Grid& grid,
        const Geometry& geometry,
        double Reynolds,
        const Flux& q_potential
    );

    /// \brief Constructor, assuming potential flow part is zero
    NavierStokesModel(
        const Grid& grid,
        const Geometry& geometry,
        double Reynolds
    );

    virtual ~NavierStokesModel();

    /// Perform initial calculations needed to use model
    void init();

    /// \brief Return true if the geometry has moving bodies
    bool isTimeDependent() const;

    /// \brief Return the number of points in the geometry
    int getNumPoints() const;

    /// \brief Return the right-hand side b of the constraint equations.
    /// Here, this is the velocity of the bodies minus the base flow velocity
    BoundaryVector getConstraints() const;
    
    /// \brief Update operators, for time-dependent models
    void updateOperators( double time );

//    /// Return a pointer to the associated Geometry
//    inline const Geometry& getGeometry() const { return _geometry; }
//
//    /// Return a pointer to the associated Grid
//    inline const Grid& getGrid() const { return _grid; }

    /// Compute gamma = B(f) as in (14)
    void B(const BoundaryVector& f, Scalar& gamma ) const;
        
    /// Compute f = C(gamma) as in (14)
    void C(const Scalar& gamma, BoundaryVector& f) const;
        
    /*! \brief Compute nonlinear terms y = N(x)
    Pure virtual function: must be overridden by subclasses.
    */
    virtual Scalar N(const State& x) const = 0;
    
    /// \brief Return the constant alpha = 1/ReynoldsNumber
    double getAlpha() const;
    
    /// Compute flux q from circulation gamma, including base flow q0
    void computeFlux(const Scalar& gamma, Flux& q ) const;

    /// \brief Compute flux q from the circulation gamma, including base flow
    void refreshState( State& x ) const;

private:    
    /*! \brief Given the circulation gamma, return the streamfunction psi.

    Assumes psi = 0 on the boundary, and does not add in potential flow solution
    */
    Scalar gammaToStreamfunction(const Scalar& gamma) const;
    BoundaryVector getBaseFlowBoundaryVelocities() const;
    void computeFluxWithoutBaseFlow(const Scalar& gamma, Flux& q ) const;

    // data
    const Grid& _grid;
    const Geometry& _geometry;
    Regularizer _regularizer;
    Flux _baseFlow;
    double _ReynoldsNumber;
    PoissonSolver _poisson;
    bool _hasBeenInitialized;
};

//! Full nonlinear Navier-Stokes equations.
class NonlinearNavierStokes : public NavierStokesModel {
public:
    NonlinearNavierStokes(
        const Grid& grid,
        const Geometry& geometry,
        double Reynolds,
        const Flux& q_potential
        ) :
        NavierStokesModel( grid, geometry, Reynolds, q_potential )
    {}

    NonlinearNavierStokes(
        const Grid& grid,
        const Geometry& geometry,
        double Reynolds
        ) :
        NavierStokesModel( grid, geometry, Reynolds )
    {}
    
    /*! \brief Compute nonlinear terms y = N(x)
    for full nonlinear Navier-Stokes equations
    */
    Scalar N(const State& x) const;
    
};

//! Navier-Stokes equations linearized about an equilibrium point.
class LinearizedNavierStokes : public NavierStokesModel {
public:
    LinearizedNavierStokes(
        const Grid& grid,
        const Geometry& geometry,
        double Reynolds,
        const State& baseFlow
        ) :
        NavierStokesModel( grid, geometry, Reynolds ),
        _x0( baseFlow )
    {}

    /*! \brief Compute nonlinear terms y = N(x)
    for linearized Navier-Stokes equations
    */
    Scalar N(const State& x) const;

private:
    State _x0;
};


//! Adjoint Navier-Stokes equations, linearized about an equilibrium point.
// For now, just make AdjointNavierStokes identical to LinearizedNavierStokes
typedef LinearizedNavierStokes AdjointNavierStokes;


} // namespace ibpm

#endif /* _NAVIERSTOKES_H_ */
