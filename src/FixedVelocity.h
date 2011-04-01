#ifndef _FIXEDVELOCITY_H_
#define _FIXEDVELOCITY_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

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
        _intX = 0.;
        _intY = 0.;
        _intTheta = 0.;
        _oldtime = 0.;
    }
    
    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, 0, theta(t), 0, 0, thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
        double dtt = time-_oldtime;
        if (dtt > 1) dtt = 0.;
        _intX = _intX + _xdot*dtt;
        _intY = _intY + _ydot*dtt;
        _intTheta = _intTheta + _thetadot*dtt;
        _oldtime = time;
        return TangentSE2( _intX, _intY, _intTheta, _xdot, _ydot, _thetadot );
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
    mutable double _intX;
    mutable double _intY;
    mutable double _intTheta;
    mutable double _oldtime;
};

} // namespace ibpm

#endif /* _FIXEDVELOCITY_H_ */
