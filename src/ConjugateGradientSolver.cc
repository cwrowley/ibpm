// ConjugateGradientSolver.cc
//
// Description:
// Implementation of the ConjugateGradientSolver class
//
// Author(s):
// Clancy Rowley
//
// Date: 8 Aug 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "NavierStokesModel.h"
#include "ProjectionSolver.h"
#include "ConjugateGradientSolver.h"

namespace ibpm {

const int MAX_ITERATIONS = 3000;

ConjugateGradientSolver::ConjugateGradientSolver(
    const Grid& grid,
    const NavierStokesModel& model,
    double beta,
    double tolerance
    ) :
    ProjectionSolver( grid, model, beta ),
    _toleranceSquared(tolerance*tolerance)
{}

// Implementation of conjugate gradient method
void ConjugateGradientSolver::Minv(
    const BoundaryVector& b,
    BoundaryVector& x
    ) {
    double alpha;
    double beta;
    double delta;
    double delta_old;
    int numIterations = 0;

    // r = b - M(x)
    BoundaryVector r = b;
    BoundaryVector q = M(x);
    r -= q;
    BoundaryVector d = r;
    delta = InnerProduct(r,r);
    
    // while error is greater than tolerance
    while ( (delta > _toleranceSquared) &&
            (numIterations < MAX_ITERATIONS) ) {
        // alpha = r^2 / <d, Md>
        q = M(d);
        alpha = delta / InnerProduct( d, q );
        x += alpha * d;
        r -= alpha * q;
        delta_old = delta;
        delta = InnerProduct( r, r );
        beta = delta / delta_old;
        d = r + beta * d;
    }
}

} // namespace ibpm
