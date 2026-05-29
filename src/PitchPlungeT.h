#ifndef _PITCHPLUNGET_H_
#define _PITCHPLUNGET_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

/*!
    \file PitchPlunge.h
    \class PitchPlunge

    \brief Subclass of Motion, for a pitching and plunging body
    \replace frequency f with period T

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 29 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class PitchPlungeT : public Motion {
public:

    /// \brief Define a Motion corresponding to sinusoidal pitching and
    /// plunging, centered about the origin:
    ///    y(t)      = A_y      sin( 2*pi * f_y t + phi_y)
    ///    \theta(t) = A_\theta sin( 2*pi * f_theta t + phi_theta)
    PitchPlungeT(
        double pitchAmplitude,
        double pitchPeriod,
        double pitchPhase,
        double plungeAmplitude,
        double plungePeriod,
        double plungePhase
        ) :
        _pitchAmp(pitchAmplitude),
        _pitchT(pitchPeriod),
        _pitchPhase(pitchPhase),
        _plungeAmp(plungeAmplitude),
        _plungeT(plungePeriod),
        _plungePhase(plungePhase) {

        double twopi = 8. * atan(1.);
        _pitchT /= twopi;
        _plungeT /= twopi;
    }

    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, y(t), theta(t), 0, ydot(t), thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
        double y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
        double ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase );
        // double theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
        // double thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );

        double theta;
        double thetadot;

        if (_pitchT == 0){
          theta    = _pitchAmp ;
          thetadot = 0;
        }
        else{
          theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
          thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );
        }
		//cout << "ydot_true..." << ydot << endl;
        return TangentSE2( 0, y, theta, 0, ydot, thetadot );
		//                 x           x_dot   x_dot = 1 ?

    }
	///////////////////////////////////
	inline double getydot(double time) const {
		double ydot = _plungeAmp / _plungeT * cos(time / _plungeT + _plungePhase);
		return ydot;
	}

  inline double gety(double time) const {
      double y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
    return y;
  }

	inline double gettheta(double time) const {
		// double theta = _pitchAmp * sin(time / _pitchT + _pitchPhase);
    double theta;

    if (_pitchT == 0){
      theta = _pitchAmp;
    }
    else{
      theta = _pitchAmp * sin(time / _pitchT + _pitchPhase);
    }
		return theta;
	}

	inline double getthetadot(double time) const {
		// double thetadot = _pitchAmp / _pitchT * cos(time / _pitchT + _pitchPhase);

    double thetadot;

    if (_pitchT == 0){
      thetadot = 0;
    }
    else{
      thetadot = _pitchAmp / _pitchT * cos(time / _pitchT + _pitchPhase);
    }
		return thetadot;
	}

	/////////////////////////////////////
    inline Motion* clone() const {
        double twopi = 8. * atan(1.);
        return new PitchPlungeT(
            _pitchAmp,
            _pitchT * twopi,
            _pitchPhase,
            _plungeAmp,
            _plungeT * twopi,
            _plungePhase
        );
    };

private:
    double _pitchAmp;
    double _pitchT;
    double _pitchPhase;
    double _plungeAmp;
    double _plungeT;
    double _plungePhase;
};

} // namespace ibpm

#endif /* _PITCHPLUNGET_H_ */
