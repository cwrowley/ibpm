#ifndef _FIXEDVELOCITY_H_
#define _FIXEDVELOCITY_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

double TTTTfv;
double INTX;
double INTY;
double INTTHETA;

/*!
    \file FixedVelocity.h
    \class FixedVelocity

    \brief Subclass of Motion, for a constant vertical velocity

    \author Steven Brunton
    \author $LastChangedBy: sbrunton $
    \date 26 Apr 2010
    \date $LastChangedDate: 2008-09-03 02:41:03 -0400 (Wed, 03 Sep 2008) $
    \version $Revision: 132 $
*/

class FixedVelocity : public Motion {
public:
    
    FixedVelocity(
        double xdot,
        double ydot,
        double thetadot
        ) :
        _xdot(xdot),
        _ydot(ydot),
        _thetadot(thetadot) {
        INTX = 0.;
        INTY = 0.;
        INTTHETA = 0.;
        TTTTfv = 0.;
    }
    
    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, 0, theta(t), 0, 0, thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
        double dtt = time-TTTTfv;
        if (dtt > 1) dtt = 0.;
        INTX = INTX + _xdot*dtt;
        INTY = INTY + _ydot*dtt;
        INTTHETA = INTTHETA + _thetadot*dtt;
        TTTTfv = time;
        cout << " INTX = " << INTX << endl;
        cout << " INTY = " << INTY << endl;
        cout << " INTTHETA = " << INTTHETA << endl;
        cout << " TTTT = " << TTTTfv << endl;
        cout << " dt = " << dtt << endl;
        return TangentSE2( INTX, INTY, INTTHETA, _xdot, _ydot, _thetadot );
    }

    inline Motion* clone() const {
        return new FixedVelocity(
            _xdot,
            _ydot,
            _thetadot
        );
    };

private:
    double _xdot;
    double _ydot;
    double _thetadot;
};

} // namespace ibpm

#endif /* _FIXEDVELOCITY_H_ */
