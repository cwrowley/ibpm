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
    
    Uses Crank-Nicolson for linear terms.  Uses the scheme given by Peyret, p. 149[3]:
    \f{align}
    Q_1&=hN(x^n)\\
    (1-\frac{h}{6}L)x_1+\frac{h}{3}Bf_1&=(1+\frac{h}{6}L)x^n+\frac{h}{3}Q_1\\
    Cx_1&=b_{n+1/3}\\
    Q_2&=-\frac{5}{9}Q_1+hN(x_1)\\
    (1-\frac{5h}{24}L)x_2+\frac{5h}{12}Bf_2&=(1+\frac{5h}{24}L)x_1+\frac{15}{16}Q_2\\
    Cx_2&=b_{n+3/4}\\
    Q_3&=-\frac{153}{128}Q_2+hN(x_2)\\
    (1-\frac{h}{8}L)x^{n+1}+\frac{h}{4}Bf^{n+1}&=(1+\frac{h}{8}L)x_2+\frac{8}{15}Q_3\\
    Cx^{n+1}&=b_{n+1} 
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
    void createAllSolvers();
    ProjectionSolver* _solver1;
    ProjectionSolver* _solver2;
    ProjectionSolver* _solver3;
    Scalar _linearTermEigenvalues1;
    Scalar _linearTermEigenvalues2;
    Scalar _linearTermEigenvalues3;
    State _x1;
    State _x2;
    Scalar _Q1;
    Scalar _Q2;
    Scalar _Q3;   
    Scalar _a;
    BoundaryVector _b;
    BoundaryVector _b0;
    double _A1, _A2, _A3; 
    double _B1, _B2, _B3;
    double _Bp1, _Bp2, _Bp3;
};

} // namespace ibpm

#endif /* _RUNGEKUTTA3_H_ */

