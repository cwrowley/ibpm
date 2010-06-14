#ifndef _LAGSTEP1_H_
#define _LAGSTEP1_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

/*!
    \file LagStep1.h
    \class LagStep1

    \brief Subclass of Motion, for a First-order lag filtered square pulse in dtheta (with area under pulse equal to AMP)

    \author Steven Brunton
    \author $LastChangedBy: sbrunton $
    \date 26 Apr 2010
    \date $LastChangedDate: 2008-09-03 02:41:03 -0400 (Wed, 03 Sep 2008) $
    \version $Revision: 132 $
*/

class LagStep1 : public Motion {
public:
    
    /// \brief Define a Motion corresponding to a first-order lag filtered square pulse in dtheta
    LagStep1(
        double AMP,
        double PW,
        double TAU,
        double t0
        ) :
        _AMP(AMP),
        _PW(PW),
        _TAU(TAU),
        _t0(t0) {
    }
    
    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, 0, theta(t), 0, 0, thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
        double pi = 4. * atan(1.);
        double AMPrad = _AMP*(pi/180); 
        double tdiff = time-_t0;
        double theta, thetadot;
        if(tdiff <= 0) {
            theta = 0.;
            thetadot = 0.;
        }
        else if( (tdiff > 0) && (tdiff <= _PW) ) {
            theta = tdiff + _TAU*exp(-tdiff/_TAU) - _TAU;
            thetadot = 1-exp(-tdiff/_TAU);
        }
        else if(tdiff > _PW) {
            theta = _PW + _TAU*exp(-tdiff/_TAU) - _TAU*exp(-(tdiff-_PW)/_TAU);
            thetadot = -1.*exp(-tdiff/_TAU) + exp(-(tdiff-_PW)/_TAU);
        }
        else {
            cerr << "LagStep1 (ERROR):  time is not compatible\n";
            exit(1);
        }
        return TangentSE2( 0, 0, theta*AMPrad/_PW, 0, 0, thetadot*AMPrad/_PW );
    }

    inline Motion* clone() const {
        return new LagStep1(
            _AMP,
            _PW,
            _TAU,
            _t0
        );
    };

private:
    double _AMP;
    double _PW;
    double _TAU;
    double _t0;
};

} // namespace ibpm

#endif /* _LAGSTEP1_H_ */
