#ifndef _ELLIPTICSOLVERMULTI_H_
#define _ELLIPTICSOLVERMULTI_H_

#include "Grid.h"
#include "Scalar.h"
#include "EllipticSolver2d.h"
#include <vector>

using std::vector;

namespace ibpm {

/*!
    \file EllipticSolver.h
    \class EllipticSolver

    \brief Solve Poisson and Helmholtz equations on a multi-domain grid (Scalar)

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 30 Sep 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class EllipticSolver {
public:
    EllipticSolver( const Grid& grid );
    virtual ~EllipticSolver() = 0;
    
    /// \brief Solve L u = f, assuming zero boundary conditions on u
    void solve( const Scalar& f, Scalar& u ) const;
    /// \brief Convenience form
    Scalar solve( const Scalar& f ) const;
protected:
    virtual EllipticSolver2d* create2dSolver( double dx ) = 0;
    void init();
private:
    int _ngrid;
    double _dx;
    vector<EllipticSolver2d *> _solvers;
};

/******************************************************************************/

/*! \class PoissonSolver
 
 \brief Solve a Poisson equation on a multi-domain grid.
 Solves L u = f, with zero boundary conditions on the outer domain of u,
 where L is the Laplacian.
 */
class PoissonSolver : public EllipticSolver {
public:
    PoissonSolver( const Grid& grid );
    EllipticSolver2d* create2dSolver( double dx );
private:
    int _nx;
    int _ny;
};

/******************************************************************************/

/*! \class HelmholtzSolver
 \brief Solve a Helmholtz equation on a multi-domain grid
 Solves (1 + alpha * L) u = f, with zero boundary conditions on u,
 where L is the Laplacian.
 */
class HelmholtzSolver : public EllipticSolver {
public:
    HelmholtzSolver( const Grid& grid, double alpha );
    EllipticSolver2d* create2dSolver( double dx );
private:
    int _nx;
    int _ny;
    double _alpha;
};
    
} // namespace ibpm

#endif /* _ELLIPTICSOLVERMULTI_H_ */
