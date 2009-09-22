#ifndef _CONJUGATEGRADIENTSOLVER_H_
#define _CONJUGATEGRADIENTSOLVER_H_

#include <math.h>

namespace ibpm {

class ProjectionSolver;
class NavierStokesModel;
class BoundaryVector;

/*!
    \file ConjugateGradientSolver.h
    \class ConjugateGradientSolver

    \brief Subclass of ProjectionSolver, in which the system Mf = b is solved iteratively, using a conjugate gradient method
    
    This class assumes the matrix M is symmetric, and iterates until a
    specified tolerance has been reached.

    \author Clancy Rowley
    \author $LastChangedBy$
    \date  9 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class ConjugateGradientSolver : public ProjectionSolver {
public:
    
    /// Constructor.  Store the tolerance as private data
    ConjugateGradientSolver(
        const Grid& grid,
        const NavierStokesModel& model,
        double beta,
        double tolerance
    );

    inline void setTolerance(double tolerance) {
        _toleranceSquared = tolerance * tolerance;
    }
    
    inline double getTolerance() {
        return sqrt(_toleranceSquared);
    }

protected:
    /// \brief Solve Mf = b for f iteratively, using a conjugate-gradient method.
    /// Assumes M is symmetric
    void Minv(
        const BoundaryVector& b,
        BoundaryVector& x
    );

private:
    double _toleranceSquared;
};

} // namespace ibpm

#endif /* _CONJUGATEGRADIENTSOLVER_H_ */

