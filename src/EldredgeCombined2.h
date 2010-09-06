#ifndef _ELDREDGECOMBINED2_H_
#define _ELDREDGECOMBINED2_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {


/*!
    \file EldredgeCombined2.h
    \class EldredgeCombined2

    \brief Subclass of Motion, for the canonical maneuver described by Eldredge

    \author Steven Brunton
    \author $LastChangedBy: sbrunton $
    \date 26 Apr 2010
    \date $LastChangedDate: 2008-09-03 02:41:03 -0400 (Wed, 03 Sep 2008) $
    \version $Revision: 132 $
*/

class EldredgeCombined2 : public Motion {
public:
    
    EldredgeCombined2(
        double AMPa,
        double a,
        double a1,
        double a2,
        double a3,
        double a4,
        double AMPb,
        double b,
        double b1,
        double b2,
        double b3,
        double b4
        ) :
        _AMPa(AMPa),
        _a(a),
        _a1(a1),
        _a2(a2),
        _a3(a3),
        _a4(a4),
        _AMPb(AMPb),
        _b(b),
        _b1(b1),
        _b2(b2),
        _b3(b3),
        _b4(b4) {
        // a variables are for plunging, b variables are for pitching
        double tta = (_a2+_a3)/2.;
        double ttb = (_b2+_b3)/2.;
        _intG = 0.;
        _oldtime = 0.;
        _maxGa = log((cosh(_a*(tta-_a1))*cosh(_a*(tta-_a4)))/(cosh(_a*(tta-_a2))*cosh(_a*(tta-_a3))));
        _maxGb = log((cosh(_b*(ttb-_b1))*cosh(_a*(ttb-_b4)))/(cosh(_b*(ttb-_b2))*cosh(_b*(ttb-_b3))));
    }
    
    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, h(t), theta(t), 0, hdot(t), thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
        double pi = 4. * atan(1.);
        double AMPrad = _AMPb*(pi/180);
        double a1d = time-_a1;
        double a2d = time-_a2;
        double a3d = time-_a3;
        double a4d = time-_a4;
        double b1d = time-_b1;
        double b2d = time-_b2;
        double b3d = time-_b3;
        double b4d = time-_b4;
        double Ga,Gb,dGdtb;
        double Garga = ((cosh(_a*a1d)*cosh(_a*a4d))/(cosh(_a*a2d)*cosh(_a*a3d)));
        double Gargb = ((cosh(_b*b1d)*cosh(_b*b4d))/(cosh(_b*b2d)*cosh(_b*b3d)));
        Ga = log(Garga);
        Gb = log(Gargb);
        dGdtb = _b*(tanh(_b*b1d)-tanh(_b*b2d)-tanh(_b*b3d)+tanh(_b*b4d));

        // The following is a hack to compute integral of G so that I can step vertical velocity...
        double dtt = time-_oldtime;
        _intG = _intG + Ga*dtt;
        _oldtime = time;

        return TangentSE2( 0, _intG*_AMPa/_maxGa, Gb*AMPrad/_maxGb, 0, Ga*_AMPa/_maxGa, dGdtb*AMPrad/_maxGb );
    }

    inline Motion* clone() const {
        return new EldredgeCombined2(
            _AMPa,
            _a,
            _a1,
            _a2,
            _a3,
            _a4,
            _AMPb,
            _b,
            _b1,
            _b2,
            _b3,
            _b4
        );
    };

private:
    double _AMPa;
    double _a;
    double _a1;
    double _a2;
    double _a3;
    double _a4;
    double _AMPb;
    double _b;
    double _b1;
    double _b2;
    double _b3;
    double _b4;
    double _maxGa;
    double _maxGb;
    mutable double _intG;
    mutable double _oldtime;
};

} // namespace ibpm

#endif /* _ELDREDGECOMBINED2_H_ */
