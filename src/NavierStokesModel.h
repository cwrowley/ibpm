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
    /*! \brief Constructor, given a specified grid, geometry, and Reynolds number
    */
    NavierStokesModel(
        const Grid& grid,
        const Geometry& geometry,
        double Reynolds,
        const Flux& q_potential
    );
    
    virtual ~NavierStokesModel();

    /// Return a pointer to the associated Geometry
    inline const Geometry* getGeometry() const { return &_geometry; }

    /// Return a pointer to the associated Grid
    inline const Grid* getGrid() const { return &_grid; }

    /// Return a pointer to the eigenvalues of the linear term L
    inline const Scalar* getLambda() const { return &_linearTermEigenvalues; }
    
    /// Transform to eigenvectors of L (discrete sin transform)
    inline Scalar S(const Scalar& g) const {
        Scalar ghat = SinTransform( g );
        return ghat;
    }

    /*! \brief Inverse transform of S.
    Note that for the discrete sin transform, the S^{-1} should equal S, but
    in some implementations (e.g. FFTW), these may differ by a normalization
    constant.
    */
    Scalar Sinv(const Scalar& ghat) const {
        Scalar g = SinTransform( ghat );
        // multiply by normalization factor
        // TODO: keep track of this normalization factor within a class
        //       containing SinTransform -- perhaps NavierStokesModel??
        int nx = _grid.getNx();
        int ny = _grid.getNy();
        double normalization = 1.0 / ( nx * ny * 4 );
        g *= normalization;
        return g;
    }
    
    /// Convenience form
    inline Scalar B(const BoundaryVector& f) const {
        Scalar gamma( _grid );
        B( f, gamma );
        return gamma;
    }

    /// Compute gamma = B(f) as in (14)
    inline void B(const BoundaryVector& f, Scalar& gamma ) const {
        Flux q = _regularizer.toGrid( f );
        gamma = Curl( q );
    }
    
    /// Compute f = C(gamma) as in (14)
    inline BoundaryVector C(const Scalar& gamma) const {
        BoundaryVector f( _geometry.getNumPoints() );
        C( gamma, f );
        return f;
    }
    
    /// Compute f = C(gamma) as in (14)
    inline void C(const Scalar& gamma, BoundaryVector& f) const {
        Flux q(_grid);
        computeFluxWithoutBaseFlow( gamma, q );
        f = _regularizer.toBoundary( q );
    }
    
    /*! \brief Compute nonlinear terms y = N(x)
    Pure virtual function: must be overridden by subclasses.
    */
    virtual Scalar nonlinear(const State& x) const = 0;
    
    /// Compute flux q from circulation gamma
    inline void computeFluxWithoutBaseFlow(const Scalar& gamma, Flux& q ) const {
        Scalar streamfunction = gammaToStreamfunction( gamma );
        q = Curl( streamfunction );
    }

    /// Compute flux q from circulation gamma, including base flow q0
    inline void computeFlux(const Scalar& gamma, Flux& q ) const {
        computeFluxWithoutBaseFlow( gamma, q );
        q += _baseFlow;
    }

    inline BoundaryVector getBaseFlowBoundaryVelocities() const {
        BoundaryVector velocity = _regularizer.toBoundary( _baseFlow );
        return velocity;
    }


protected:
    
    /*! \brief Given the circulation gamma, return the streamfunction psi.

    TODO: Assumptions about boundary conditions??
    */
    inline Scalar gammaToStreamfunction(const Scalar& gamma) const {
        Scalar psi = S( gamma );
        psi *= _eigGammaToStreamfunction;
        psi = Sinv( psi );
        return psi;
    }
    
private:
    const Grid& _grid;
    const Geometry& _geometry;
    Regularizer _regularizer;
    Scalar _linearTermEigenvalues;
    Scalar _eigGammaToStreamfunction;
    Flux _baseFlow;
    double _ReynoldsNumber;
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
        NavierStokesModel( grid, geometry, Reynolds, q_potential ) {}


    /*! \brief Compute nonlinear terms y = N(x)
    for full nonlinear Navier-Stokes equations
    */
    inline Scalar nonlinear(const State& x) const {
        Flux v = CrossProduct( x.q, x.gamma );
        Scalar g = Curl( v );
        return g;
    };
    
};

//! Navier-Stokes equations linearized about an equilibrium point.
class LinearizedNavierStokes : public NavierStokesModel {};

//! Adjoint Navier-Stokes equations, linearized about an equilibrium point.
class AdjointNavierStokes : public NavierStokesModel {};

} // namespace ibpm

#endif /* _NAVIERSTOKES_H_ */
