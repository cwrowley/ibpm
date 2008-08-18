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
	    ) :
	    _grid(grid),
	    _geometry(geometry),
        _regularizer( grid, geometry ),
        _linearTermEigenvalues( grid ),
        _inverseLaplacianEigenvalues( grid ),
        _baseFlow(q_potential),
	    _ReynoldsNumber(Reynolds)
	    {
        
        // calculate eigenvalues of Laplacian
        int nx = grid.getNx();
        int ny = grid.getNy();    
        const double pi = 4. * atan(1.);
        Scalar eigLaplacian(grid);
        eigLaplacian = 1.;
        // Loop over only interior points
        for (int i=1; i < nx; ++i ) {
            for (int j=1; j < ny; ++j ) {
                eigLaplacian(i,j) = 2. * ( cos( (pi * i) / nx ) +
                                           cos( (pi * j) / ny ) - 2. ) / 
                                    ( _grid.getDx() * _grid.getDx() );
            }
        }
            
        _inverseLaplacianEigenvalues = 1. / eigLaplacian;
        
        // calculate linear term
        double beta = 1. / Reynolds;
        _linearTermEigenvalues = beta * eigLaplacian;

        // Update regularizer
        _regularizer.update();
    }

	~NavierStokesModel() {}

	/// Return a pointer to the associated Geometry
	const Geometry* getGeometry() const { return &_geometry; }

	/// Return a pointer to the associated Grid
	const Grid* getGrid() const { return &_grid; }

	/// Return a pointer to the eigenvalues of the linear term L
	const Scalar* getLambda() const { return &_linearTermEigenvalues; }
	
	/// Transform to eigenvectors of L (discrete sin transform)
	inline Scalar S(const Scalar& g) const {
        Scalar ghat = sinTransform( g );
        return ghat;
	}

	/*! \brief Inverse transform of S.
	Note that for the discrete sin transform, the S^{-1} should equal S, but
	in some implementations (e.g. FFTW), these may differ by a normalization
	constant.
	*/
	Scalar Sinv(const Scalar& ghat) const {
        Scalar g = sinTransform( ghat );
        // multiply by normalization factor
        // TODO: keep track of this normalization factor within a class
        //       containing sinTransform -- perhaps NavierStokesModel??
        int nx = _grid.getNx();
        int ny = _grid.getNy();
        double normalization = 1.0 / ( nx * ny * 4 );
        g *= normalization;
        return g;
	}
	
	/// Compute gamma = B(f) as in (14)
	inline Scalar B(const BoundaryVector& f) const {
        Flux q = _regularizer.toGrid( f );
        Scalar gamma = curl( q );
        return gamma;
	}
	
	/// Compute f = C(gamma) as in (14)
	inline BoundaryVector C(const Scalar& gamma) const {
        Flux q(_grid);
        computeFlux( gamma, q );
        BoundaryVector f = _regularizer.toBoundary( q );
        return f;
	}
	
	/*! \brief Compute nonlinear terms y = N(x)
	Pure virtual function: must be overridden by subclasses.
	*/
	virtual Scalar nonlinear(const State& x) const = 0;
	
	/// Compute flux q from circulation gamma
    void computeFlux( const Scalar& gamma, Flux& q ) const {
        Scalar streamfunction = inverseLaplacian( gamma );
        q = curl( streamfunction );
        q += _baseFlow;
    }

protected:
    
    /*! \brief Return the inverse Laplacian of gamma
    Assumptions about boundary conditions??
    */
    inline Scalar inverseLaplacian(const Scalar& gamma) const {
        Scalar ghat = S( gamma );
        ghat *= _inverseLaplacianEigenvalues;
        ghat = Sinv( ghat );
        return ghat;
    }
    
	/*! \brief Compute bilinear term, used by subclasses
	(e.g. y = N(q) = bilinear(q,q) )
	*/
	Scalar bilinear(const State& x1, const State& x2);
	
private:
	const Geometry& _geometry;
	const Grid& _grid;
    Regularizer _regularizer;
	Scalar _linearTermEigenvalues;
    Scalar _inverseLaplacianEigenvalues;
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
        Flux v = crossproduct( x.q, x.gamma );
        Scalar g = curl( v );
        return g;
	};
    
};

//! Navier-Stokes equations linearized about an equilibrium point.
class LinearizedNavierStokes : public NavierStokesModel {};

//! Adjoint Navier-Stokes equations, linearized about an equilibrium point.
class AdjointNavierStokes : public NavierStokesModel {};


#endif /* _NAVIERSTOKES_H_ */
