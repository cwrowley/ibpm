#ifndef _LAGSTEP2_H_
#define _LAGSTEP2_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

/*!
    \file LagStep2.h
    \class LagStep2

    \brief Subclass of Motion, for a second-order lag filtered square pulse in dtheta (with area under pulse equal to AMP)

    \author Steven Brunton
    \author $LastChangedBy: sbrunton $
    \date 26 Apr 2010
    \date $LastChangedDate: 2008-09-03 02:41:03 -0400 (Wed, 03 Sep 2008) $
    \version $Revision: 132 $
*/

class LagStep2 : public Motion {
public:
    
    /// \brief Define a Motion corresponding to a second-order lag filtered square pulse in dtheta
    LagStep2(
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
            theta = tdiff + (2*_TAU+tdiff)*exp(-tdiff/_TAU) - 2*_TAU; 
            thetadot = 1.-(1.+tdiff/_TAU)*exp(-tdiff/_TAU);
        }
        else if(tdiff > _PW) {
            theta = _PW + (2*_TAU+tdiff)*exp(-tdiff/_TAU) + (_PW-2*_TAU)*exp(-(tdiff-_PW)/_TAU) - tdiff*exp(-(tdiff-_PW)/_TAU); 
            thetadot = -(1.+tdiff/_TAU)*exp(-tdiff/_TAU) + (1.+(tdiff-_PW)/_TAU)*exp(-(tdiff-_PW)/_TAU);
        }
        else {
            cerr << "LagStep2 (ERROR):  time is not compatible\n";
            exit(1);
        }
        return TangentSE2( 0, 0, theta*AMPrad/_PW, 0, 0, thetadot*AMPrad/_PW );
    }

    inline Motion* clone() const {
        return new LagStep2(
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

#endif /* _LAGSTEP2_H_ */
