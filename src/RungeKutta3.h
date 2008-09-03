#ifndef _RUNGEKUTTA3_H_
#define _RUNGEKUTTA3_H_

#include "TimeStepper.h"
#include "State.h"
#include <string>
using std::string;

namespace ibpm {

/*!
    \file RungeKutta3.h
    \class RungeKutta3

    \brief Timestepper using RK3 for nonlinear terms
    
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

class RungeKutta3 : public TimeStepper {
public:

    /// Instantiate an RK3 solver
    RungeKutta3( NavierStokesModel& model, double timestep );

    /// Destructor
    ~RungeKutta3();

    void init();
    bool load(const std::string& basename);
    bool save(const std::string& basename);

    /// Advance the state forward one step, using RK3
    void advance(State& x);
private:
    ProjectionSolver* _solver1;
    ProjectionSolver* _solver2;
    ProjectionSolver* _solver3;
    Scalar _linearTermEigenvalues1;
    Scalar _linearTermEigenvalues2;
    Scalar _linearTermEigenvalues3;
    State _x1;
    State _x2;
};

} // namespace ibpm

#endif /* _RUNGEKUTTA3_H_ */

