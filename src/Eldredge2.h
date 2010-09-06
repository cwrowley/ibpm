#ifndef _ELDREDGE2_H_
#define _ELDREDGE2_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {


/*!
    \file Eldredge2.h
    \class Eldredge2

    \brief Subclass of Motion, for the canonical maneuver described by Eldredge

    \author Steven Brunton
    \author $LastChangedBy: sbrunton $
    \date 26 Apr 2010
    \date $LastChangedDate: 2008-09-03 02:41:03 -0400 (Wed, 03 Sep 2008) $
    \version $Revision: 132 $
*/

class Eldredge2 : public Motion {
public:
    
    Eldredge2(
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
        cout << " got her"<< endl;
        _intG = 0.;
        _oldtime = 0.;
        cout << " not here" << endl;
        _maxG = log((cosh(_a*(tt-_t1))*cosh(_a*(tt-_t4)))/(cosh(_a*(tt-_t2))*cosh(_a*(tt-_t3))));
//        _maxG = log((cosh(_a*(tt-t1)))/(cosh(_a*(tt-_t2))));
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
        //G = log((cosh(_a*t1d)*cosh(_a*t4d))/(cosh(_a*t2d)*cosh(_a*t3d)));
        dGdt = _a*(tanh(_a*t1d)-tanh(_a*t2d)-tanh(_a*t3d)+tanh(_a*t4d));

        // The following is a hack to compute integral of G so that I can step vertical velocity...
        double dtt = time-_oldtime;
        _intG = _intG + G*dtt;
        _oldtime = time;
        cout << "dtt = " << dtt << endl;
        cout << "intG = " << _intG << endl;
        cout << "G = " << G*_AMP/_maxG << endl;
        cout << "intG*stuff = " << _intG*_AMP/_maxG << endl;
        return TangentSE2( 0, _intG*_AMP/_maxG, 0, 0, G*_AMP/_maxG, 0 );
    }

    inline Motion* clone() const {
        return new Eldredge2(
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
    mutable double _intG;
    mutable double _oldtime;
};

} // namespace ibpm

#endif /* _ELDREDGE2_H_ */
