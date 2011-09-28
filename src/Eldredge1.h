#ifndef _ELDREDGE1_H_
#define _ELDREDGE1_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {
    
/*!
    \file Eldredge1.h
    \class Eldredge1

    \brief Subclass of Motion, for the canonical maneuver described by Eldredge

    \author Steven Brunton
    \author $LastChangedBy: sbrunton $
    \date 26 Apr 2010
    \date $LastChangedDate: 2008-09-03 02:41:03 -0400 (Wed, 03 Sep 2008) $
    \version $Revision: 132 $
*/

class Eldredge1 : public Motion {
public:
    
    Eldredge1(
        double AMP,
        double a,
        double t1,
        double t2,
        double t3,
        double t4
        ) :
        _AMP(AMP),
        _a(a),
        _t1(t1),
        _t2(t2),
        _t3(t3),
        _t4(t4) {
        double tt = (_t2+_t3)/2.;
        _maxG = log((cosh(_a*(tt-_t1))*cosh(_a*(tt-_t4)))/(cosh(_a*(tt-_t2))*cosh(_a*(tt-_t3))));
    }
    
    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, 0, theta(t), 0, 0, thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
        double t1d = time-_t1;
        double t2d = time-_t2;
        double t3d = time-_t3;
        double t4d = time-_t4;
        double G,dGdt;
        double Garg = ((cosh(_a*t1d)*cosh(_a*t4d))/(cosh(_a*t2d)*cosh(_a*t3d)));
        G = log(Garg);
        dGdt = _a*(tanh(_a*t1d)-tanh(_a*t2d)-tanh(_a*t3d)+tanh(_a*t4d));

        return TangentSE2( 0, G*_AMP/_maxG, 0, 0, dGdt*_AMP/_maxG, 0 );
    }

    inline Motion* clone() const {
        return new Eldredge1(
            _AMP,
            _a,
            _t1,
            _t2,
            _t3,
            _t4
        );
    };

private:
    double _AMP;
    double _a;
    double _t1;
    double _t2;
    double _t3;
    double _t4;
    double _maxG;
};

} // namespace ibpm

#endif /* _ELDREDGE1_H_ */
