#ifndef _PITCHPLUNGET_ADJ_CHECK_H_
#define _PITCHPLUNGET_ADJ_CHECK_H_

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

class PitchPlungeT_adj_check : public Motion {
public:

    /// \brief Define a Motion corresponding to sinusoidal pitching and
    /// plunging, centered about the origin:
    ///    y(t)      = A_y      sin( 2*pi * f_y t + phi_y)
    ///    \theta(t) = A_\theta sin( 2*pi * f_theta t + phi_theta)
    PitchPlungeT_adj_check(
        double pitchAmplitude,
        double pitchPeriod,
        double pitchPhase,
        double plungeAmplitude,
        double plungePeriod,
        double plungePhase,
        double offsety,
        double offsetlocationy,
        double dt
        ) :
        _pitchAmp(pitchAmplitude),
        _pitchT(pitchPeriod),
        _pitchPhase(pitchPhase),
        _plungeAmp(plungeAmplitude),
        _plungeT(plungePeriod),
        _plungePhase(plungePhase),
        _offsety(offsety),
        _offsetlocationy(offsetlocationy),
        _dt(dt) {

        double twopi = 8. * atan(1.);
        _pitchT /= twopi;
        _plungeT /= twopi;
    }

    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, y(t), theta(t), 0, ydot(t), thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
      double step = time / _dt;
      double y;
      double ydot;
      if ((int)step == (int)_offsetlocationy){ // remove int!!!!!!!
        y    = _plungeAmp * sin( time / _plungeT + _plungePhase ) + _offsety;
        ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase ) - _offsety / _dt;
        cout << " case1 ..." << endl;
      }
     else if ((int)step == _offsetlocationy - 1){
        y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
        ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase ) + _offsety / _dt;
        cout << " case2 ..." << endl;
      }
      else{
        y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
        ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase );
        cout << " case3 ..." << endl;
      }

        double theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
        double thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );
		//cout << "ydot_true..." << ydot << endl;
    cout << " step... " << step << endl;
    cout << " offsetlocationy... " << (int)_offsetlocationy << endl;
    cout << " y... " << y << endl;
    cout << " ydot..." << ydot << endl;
        return TangentSE2( 0, y, theta, 0, ydot, thetadot );
		//                 x           x_dot

    }
	///////////////////////////////////
	inline double getydot(double time) const {
    double step = time / _dt;
    double ydot;
    if ((int)step == (int)_offsetlocationy){
      ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase ) - _offsety / _dt;
    }
    if ((int)step == _offsetlocationy - 1){
      ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase ) + _offsety / _dt;
    }
    else{
      ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase );
    }
		return ydot;
	}

  inline double gety(double time) const {
    double step = time / _dt;
    double y;
    if ((int)step == (int)_offsetlocationy){
      y    = _plungeAmp * sin( time / _plungeT + _plungePhase ) + _offsety;
    }
    if ((int)step == _offsetlocationy - 1){
      y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
    }
    else{
      y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
    }

    return y;
  }

	inline double gettheta(double time) const {
      double theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
		return theta;
	}

	inline double getthetadot(double time) const {
      double thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );
		return thetadot;
	}

	/////////////////////////////////////
    inline Motion* clone() const {
        double twopi = 8. * atan(1.);
        return new PitchPlungeT_adj_check(
            _pitchAmp,
            _pitchT * twopi,
            _pitchPhase,
            _plungeAmp,
            _plungeT * twopi,
            _plungePhase,
            _offsety,
            _offsetlocationy,
            _dt
        );
    };

private:
    double _pitchAmp;
    double _pitchT;
    double _pitchPhase;
    double _plungeAmp;
    double _plungeT;
    double _plungePhase;
    double _offsety;
    double _offsetlocationy;
    double _dt;
};

} // namespace ibpm

#endif /* _PITCHPLUNGET_H_ */
