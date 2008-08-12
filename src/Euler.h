#ifndef _EULER_H_
#define _EULER_H_

/*!
    \file Euler.h
    \class Euler

    \brief Timestepper using explicit Euler for nonlinear terms
    
    Uses Crank-Nicolson for linear terms.  In particular, advance() returns the solution of
    \f{align}
    (1 - \frac{h}{2}L)\gamma^{n+1} + hBf &=
        (1 + \frac{h}{2}L) \gamma^n + h N(q^n)\\
    C\gamma^{n+1} &= b_{n+1}
    \f}
    

    \author Clancy Rowley
    \author $LastChangedBy$
    \date  2 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/
#include "TimeStepper.h"

class Euler : public TimeStepper {
public:

    /// Instantiate an Euler solver
    Euler( const NavierStokesModel& model, double timestep );

    /// Advance the state forward one step, using explicit Euler
    void advance(State& x);

private:
    ProjectionSolver* _solver;
};


#endif /* _EULER_H_ */

