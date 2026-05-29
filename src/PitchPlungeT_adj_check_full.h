#ifndef _PITCHPLUNGET_ADJ_CHECK_FULL_H_
#define _PITCHPLUNGET_ADJ_CHECK_FULL_H_

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

class PitchPlungeT_adj_check_full : public Motion {
public:

    /// \brief Define a Motion corresponding to sinusoidal pitching and
    /// plunging, centered about the origin:
    ///    y(t)      = A_y      sin( 2*pi * f_y t + phi_y)
    ///    \theta(t) = A_\theta sin( 2*pi * f_theta t + phi_theta)
    PitchPlungeT_adj_check_full(
        double pitchAmplitude,
        double pitchPeriod,
        double pitchPhase,
        double plungeAmplitude,
        double plungePeriod,
        double plungePhase,
        double offsettheta,
        double offsetlocationtheta,
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
        _offsettheta(offsettheta),
        _offsetlocationtheta(offsetlocationtheta),
        _offsety(offsety),
        _offsetlocationy(offsetlocationy),
        _dt(dt) {

        double twopi = 8. * atan(1.);
        _pitchT /= twopi;
        _plungeT /= twopi;
        _offsetlocationy = round(_offsetlocationy);
        _offsetlocationtheta = round(_offsetlocationtheta);
    }

    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, y(t), theta(t), 0, ydot(t), thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
      double twopi = 8. * atan(1.);
      double step = time / _dt;
      step = std::round(step);
      double period = _plungeT * twopi / _dt;
      if (step > period){
        step = std::remainder(step,period);
      }
      if (step < 0){
        step = period + step;
      }

      double y;
      double ydot;
      double theta;
      double thetadot;
      if ((int)step == (int)_offsetlocationy){ // remove int!!!!!!!
        y    = _plungeAmp * sin( time / _plungeT + _plungePhase ) + _offsety;
        ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase ) - _offsety / _dt;
        // theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
        // thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );

        if (_pitchT == 0 ){
          theta    = _pitchAmp;
          thetadot = 0;
        }
        else{
          theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
          thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );
        }

        cout << " case1 ..." << endl;
      }
     else if ((int)step == (int)_offsetlocationy - 1){
        y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
        ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase ) + _offsety / _dt;
        // theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
        // thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );

        if (_pitchT == 0 ){
          theta    = _pitchAmp;
          thetadot = 0;
        }
        else{
          theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
          thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );
        }

        cout << " case2 ..." << endl;
      }
     else if ((int)step == (int)_offsetlocationtheta){
       y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
       ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase );
       // theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase ) + _offsettheta;
       // thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase ) - _offsettheta / _dt;

       if (_pitchT == 0 ){
         theta    = _pitchAmp + _offsettheta;
         thetadot = 0 - _offsettheta / _dt;
       }
       else{
         theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase ) + _offsettheta;
         thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase ) - _offsettheta / _dt;
       }

       cout << " case3 ..." << endl;
     }
     else if ((int)step == (int)_offsetlocationtheta - 1){
       y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
       ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase );
       // theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
       // thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase ) + _offsettheta / _dt;

       if(_pitchT == 0){
         theta    = _pitchAmp;
         thetadot = 0 + _offsettheta / _dt;
       }
       else{
         theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
         thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase ) + _offsettheta / _dt;
       }
       cout << " case4 ..." << endl;
     }

      else{
        y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
        ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase );
        // theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
        // thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );

        if (_pitchT == 0){
          theta    = _pitchAmp;
          thetadot = 0;
        }
        else{
          theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
          thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );
        }
        cout << " case5 ..." << endl;
      }


		//cout << "ydot_true..." << ydot << endl;
    // cout << " step... " << (int)step << endl;
    cout << " time... " << time << endl;
    cout << " step... " << (int)step << endl;
    cout << " offsetlocationy... " << (int)_offsetlocationy << endl;
    cout << " y... " << y << endl;
    cout << " ydot..." << ydot << endl;
    cout << " offsetlocationtheta... " << (int)_offsetlocationtheta << endl;
    cout << " theta... " << theta << endl;
    cout << " thetadot..." << thetadot << endl;
    cout << " period..." << period << endl;
        return TangentSE2( 0, y, theta, 0, ydot, thetadot );
		//                     x           x_dot

    }
	///////////////////////////////////
	inline double getydot(double time) const {
    double twopi = 8. * atan(1.);
    double step= time / _dt;
    step = std::round(step);
    double period = _plungeT * twopi / _dt;
    if (step > period){
      step = std::remainder(step,period);
    }
    if (step < 0){
      step = period + step;
    }

    double ydot;
    if ((int)step == (int)_offsetlocationy){
      ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase ) - _offsety / _dt;
    }
    else if ((int)step == (int)_offsetlocationy - 1){
      ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase ) + _offsety / _dt;
    }
    else{
      ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase );
    }
		return ydot;
	}

  inline double gety(double time) const {
    double twopi = 8. * atan(1.);
    double step= time / _dt;
    step = std::round(step);
    double period = _plungeT * twopi / _dt;
    if (step > period){
      step = std::remainder(step,period);
    }
    if (step < 0){
      step = period + step;
    }

    double y;
    if ((int)step == (int)_offsetlocationy){
      y    = _plungeAmp * sin( time / _plungeT + _plungePhase ) + _offsety;
    }
    else if ((int)step == (int)_offsetlocationy - 1){
      y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
    }
    else{
      y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
    }

    return y;
  }

	inline double gettheta(double time) const {
    double twopi = 8. * atan(1.);
    double step= time / _dt;
    step = std::round(step);
    double period = _plungeT * twopi / _dt;
    if (step > period){
      step = std::remainder(step,period);
    }
    if (step < 0){
      step = period + step;
    }

    double theta;

   if ((int)step == (int)_offsetlocationtheta){
     // theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase ) + _offsettheta;

     if(_pitchT == 0){
       theta    = _pitchAmp + _offsettheta;
     }
     else{
       theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase ) + _offsettheta;
     }
   }
   else if ((int)step == (int)_offsetlocationtheta - 1){
     // theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );

     if (_pitchT == 0){
       theta    = _pitchAmp;
     }
     else{
       theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
     }
   }

   else{
      // theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );

      if(_pitchT == 0){
        theta    = _pitchAmp;
      }
      else{
        theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
      }
    }
		return theta;
	}

	inline double getthetadot(double time) const {
    double twopi = 8. * atan(1.);
    double step= time / _dt;
    step = std::round(step);
    double period = _plungeT * twopi / _dt;
    if (step > period){
      step = std::remainder(step,period);
    }
    if (step < 0){
      step = period + step;
    }

    double thetadot;
   if ((int)step == (int)_offsetlocationtheta){
     // thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase ) - _offsettheta / _dt;

     if(_pitchT == 0){
       thetadot = 0 - _offsettheta / _dt;
     }
     else{
       thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase ) - _offsettheta / _dt;
     }
   }
   else if ((int)step == (int)_offsetlocationtheta - 1){
     // thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase ) + _offsettheta / _dt;
     if (_pitchT == 0){
       thetadot = 0 + _offsettheta / _dt;
     }
     else{
       thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase ) + _offsettheta / _dt;
     }
   }

    else{
      // thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );

      if (_pitchT == 0){
        thetadot = 0;
      }
      else{
        thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );
      }
    }
		return thetadot;
	}

	/////////////////////////////////////
    inline Motion* clone() const {
        double twopi = 8. * atan(1.);
        return new PitchPlungeT_adj_check_full(
            _pitchAmp,
            _pitchT * twopi,
            _pitchPhase,
            _plungeAmp,
            _plungeT * twopi,
            _plungePhase,
            _offsettheta,
            _offsetlocationtheta,
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
    double _offsettheta;
    double _offsetlocationtheta;
    double _offsety;
    double _offsetlocationy;
    double _dt;
};

} // namespace ibpm

#endif /* _PITCHPLUNGET_H_ */
