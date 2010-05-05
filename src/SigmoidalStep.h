#ifndef _SIGMOIDALSTEP_H_
#define _SIGMOIDALSTEP_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

/*!
    \file SigmoidalStep.h
    \class SigmoidalStep

    \brief Subclass of Motion, for a sigmoidal step in angle-of-attack

    \author Steven Brunton
    \author $LastChangedBy: sbrunton $
    \date 26 Apr 2010
    \date $LastChangedDate: 2008-09-03 02:41:03 -0400 (Wed, 03 Sep 2008) $
    \version $Revision: 132 $
*/

class SigmoidalStep : public Motion {
public:
    
    /// \brief Define a Motion corresponding to a sigmoidal pitch-up
    /// centered about the origin:
    ///    \alpha(t) = \text{AMP}\cdot \frac{1}{2} \left(1+\erf\left(12\frac{t}{\text{DUR}}-6\right)\right)
    SigmoidalStep(
        double AMP,
        double DUR,
        double t0
        ) :
        _AMP(AMP),
        _DUR(DUR),
        _t0(t0) {
    }
    
    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, 0, theta(t), 0, 0, thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
        double pi = 4. * atan(1.);
        double tdiff = time-_t0;
        double arg = ((12.*tdiff/_DUR)-6.);
        double coeff = 12./(_DUR*sqrt(pi));
        double nerf = erf(arg);
        double sig = (1./2)*_AMP*(pi/180)*(1.+nerf); 
        double sigdot = _AMP*(pi/180)*coeff*exp(-arg*arg);
        double theta    = sig;
        double thetadot = sigdot;
        return TangentSE2( 0, 0, theta, 0, 0, thetadot );
    }

    inline Motion* clone() const {
        return new SigmoidalStep(
            _AMP,
            _DUR,
            _t0
        );
    };

private:
    double _AMP;
    double _DUR;
    double _t0;
};

} // namespace ibpm

#endif /* _SIGMOIDALSTEP_H_ */
