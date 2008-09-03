#ifndef _ADAMSBASHFORTH_H_
#define _ADAMSBASHFORTH_H_

#include "TimeStepper.h"
#include "State.h"
#include "Geometry.h"
#include "Grid.h"

namespace ibpm {

/*!
    \file AdamsBashforth.h
    \class AdamsBashforth

    \brief Timestepper using Adams-Bahsforth for nonlinear terms
    
    Uses Crank-Nicolson for linear terms.  Uses the scheme given by Peyret, p. 148[3], for alpha=1, beta=1/2:  
    \f{align}
    (1 - \frac{h}{2}L)\gamma^{n+1} + hBf &=
        (1+\frac{h}{2}L)\gamma^n + \frac{h}{2}(3N(q^n)-N(q^{n-1}))\\
    C\gamma^{n+1} &= b_{n+1} 
    \f}  

    \author Steve Brunton
    \author $LastChangedBy: sbrunton $
    \date  28 Aug 2008
    \date $LastChangedDate: 2008-08-13 11:11:32 -0400 (Thu, 28 Aug 2008) $
    \version $Revision: 105 $
*/

class AdamsBashforth : public TimeStepper {
public:

    /// Instantiate an AB solver
    AdamsBashforth( NavierStokesModel& model, double timestep );

    /// Destructor
    ~AdamsBashforth();

    /// Set previous state _xold
    inline void setPreviousState(const State& x) {
        _xold = x;
        _oldSaved = true;
    }

    inline void setTempState(const State& x) {
        _xtemp = x;
    }
 
    void init();
    bool load(const std::string& basename);
    bool save(const std::string& basename);

    /// Advance the state forward one step, using AB
    void advance(State& x);
private:
    ProjectionSolver* _solver;
    Scalar _linearTermEigenvalues;
    State _xold;
    State _xtemp;  
    bool _oldSaved;
};

} // namespace ibpm

#endif /* _ADAMSBASHFORTH_H_ */

