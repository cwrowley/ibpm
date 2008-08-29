#ifndef _RUNGEKUTTA2_H_
#define _RUNGEKUTTA2_H_

#include "TimeStepper.h"
#include "State.h"
#include "Geometry.h"
#include "Grid.h"

/*!
    \file RungeKutta2.h
    \class RungeKutta2

    \brief Timestepper using RK2 for nonlinear terms
    
    Uses Crank-Nicolson for linear terms.  Uses the scheme given by Peyret, p. 148[3], for alpha=1, beta=1/2:  
    \f{align}
    (1 - \frac{h}{2}L)\gamma_1 + hBf_1 &=
        (1+\frac{h}{2}L)\gamma^n + h N(q^n)\\
    C\gamma_1 &= b_{n+1} \\
    (1 - \frac{h}{2}L)\gamma^{n+1} + hBf^{n+1} &=
        (1+\frac{h}{2}L)\gamma^n + \frac{h}{2}(N(q^n)+N(q_1))\\
    C\gamma^{n+1} &= b_{n+1}
    \f}  

    \author Steve Brunton
    \author $LastChangedBy: sbrunton $
    \date  28 Aug 2008
    \date $LastChangedDate: 2008-08-13 11:11:32 -0400 (Thu, 28 Aug 2008) $
    \version $Revision: 105 $
*/

class RungeKutta2 : public TimeStepper {
public:

    /// Instantiate an RK2 solver
    RungeKutta2( const NavierStokesModel& model, double timestep );

    /// Advance the state forward one step, using RK2
    void advance(State& x);
private:
    ProjectionSolver* _solver;
    Scalar _linearTermEigenvalues;
    State _x1;
};


#endif /* _RUNGEKUTTA2_H_ */

