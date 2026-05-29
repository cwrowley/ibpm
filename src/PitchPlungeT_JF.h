#ifndef _PITCHPLUNGET_JF_H_
#define _PITCHPLUNGET_JF_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>
#include <boost/math/special_functions.hpp>

using namespace boost;

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

class PitchPlungeT_JF : public Motion {
public:

    /// \brief Define a Motion corresponding to sinusoidal pitching and
    /// plunging, centered about the origin:
    ///    y(t)      = A_y      sin( 2*pi * f_y t + phi_y)
    ///    \theta(t) = A_\theta sin( 2*pi * f_theta t + phi_theta)
    PitchPlungeT_JF(
        double pitchAmplitude,

        double pitchPhase,
        double plungeAmplitude,

        double plungePhase,

        double dt,
        double m




        ) :
        _pitchAmp(pitchAmplitude),

        _pitchPhase(pitchPhase),
        _plungeAmp(plungeAmplitude),

        _plungePhase(plungePhase),

        _dt(dt),
        _m(m){


    }



    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, y(t), theta(t), 0, ydot(t), thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
    //////////////
      // const double PI = 4. * atan(1.);
      // double u[1000];
      // for (int i=0; i<999; i++){
      //   u[i] = double(i)/100;
      // }
      // double k=0.99;
      // double sn = math::jacobi_sn(k,0.0);
      // for (int i = 0; i<1000; i++){
      //   double sn = math::jacobi_sn(k,u[i]);
      //   double K = math::ellint_1(k);
      //   cout << "sn" << sn << endl;
      //   cout << "k" << K << endl;
      // }
    /////////////////
      double u;
      double u_plus;
      double y;
      double y_plus;
      double ydot;
      double theta;
      double theta_plus;
      double thetadot;

      if (_m >= 0){
        u = time;
        u_plus = time + _dt;
        // Jacobi function
        y = _plungeAmp * math::jacobi_sn(_m,u);
        theta = _pitchAmp * math::jacobi_sn(_m,u);
        y_plus = _plungeAmp * math::jacobi_sn(_m,u_plus);
        theta_plus = _pitchAmp * math::jacobi_sn(_m,u_plus);

        ydot = (y_plus - y) / _dt;
        thetadot = (theta_plus - theta) / _dt;
      }
      else {
        double mu = abs(_m) / (1 + abs(_m));
        double mu1 = 1 / (1 + abs(_m));
        u = time;
        u_plus = time + _dt;

        double sn_temp = math::jacobi_sn(mu,u);
        double dn = math::jacobi_dn(mu,u);
        double sn = sqrt(mu1) * sn_temp / dn;

        double sn_temp_plus = math::jacobi_sn(mu,u_plus);
        double dn_plus = math::jacobi_dn(mu,u_plus);
        double sn_plus = sqrt(mu1) * sn_temp_plus / dn_plus;

        y = _plungeAmp * sn;
        theta = _pitchAmp * sn;
        y_plus = _plungeAmp * sn_plus;
        theta_plus = _pitchAmp * sn_plus;

        ydot = (y_plus - y) / _dt;
        thetadot = (theta_plus - theta) / _dt;

      }

      //




    /////////////////
        // double y    = _plungeAmp * sin( time / _plungeT + _plungePhase );
        // double ydot = _plungeAmp / _plungeT * cos( time /_plungeT + _plungePhase );
        // double theta    = _pitchAmp * sin( time / _pitchT + _pitchPhase );
        // double thetadot = _pitchAmp / _pitchT * cos( time / _pitchT + _pitchPhase );
		//cout << "ydot_true..." << ydot << endl;
        return TangentSE2( 0, y, theta, 0, ydot, thetadot );
		//                 x           x_dot

    }
	///////////////////////////////////
	inline double getydot(double time) const {
    /////////////////
      double u;
      double u_plus;
      double y;
      double y_plus;
      double ydot;


      if (_m >= 0){
        u = time;
        u_plus = time + _dt;
        // Jacobi function
        y = _plungeAmp * math::jacobi_sn(_m,u);

        y_plus = _plungeAmp * math::jacobi_sn(_m,u_plus);

        ydot = (y_plus - y) / _dt;

      }
      else {
        double mu = abs(_m) / (1 + abs(_m));
        double mu1 = 1 / (1 + abs(_m));
        u = time;
        u_plus = time + _dt;

        double sn_temp = math::jacobi_sn(mu,u);
        double dn = math::jacobi_dn(mu,u);
        double sn = sqrt(mu1) * sn_temp / dn;

        double sn_temp_plus = math::jacobi_sn(mu,u_plus);
        double dn_plus = math::jacobi_dn(mu,u_plus);
        double sn_plus = sqrt(mu1) * sn_temp_plus / dn_plus;

        y = _plungeAmp * sn;

        y_plus = _plungeAmp * sn_plus;

        ydot = (y_plus - y) / _dt;


      }

      //
		return ydot;
	}

	inline double gettheta(double time) const {
    /////////////////
      double u;

      double theta;


      if (_m >= 0){
        u = time;
        // Jacobi function

        theta = _pitchAmp * math::jacobi_sn(_m,u);

      }
      else {
        double mu = abs(_m) / (1 + abs(_m));
        double mu1 = 1 / (1 + abs(_m));
        u = time;

        double sn_temp = math::jacobi_sn(mu,u);
        double dn = math::jacobi_dn(mu,u);
        double sn = sqrt(mu1) * sn_temp / dn;



        theta = _pitchAmp * sn;

      }

      //
		return theta;
	}

	inline double getthetadot(double time) const {
    /////////////////
      double u;
      double u_plus;

      double theta;
      double theta_plus;
      double thetadot;

      if (_m >= 0){
        u = time;
        u_plus = time + _dt;
        // Jacobi function

        theta = _pitchAmp * math::jacobi_sn(_m,u);

        theta_plus = _pitchAmp * math::jacobi_sn(_m,u_plus);


        thetadot = (theta_plus - theta) / _dt;
      }
      else {
        double mu = abs(_m) / (1 + abs(_m));
        double mu1 = 1 / (1 + abs(_m));
        u = time;
        u_plus = time + _dt;

        double sn_temp = math::jacobi_sn(mu,u);
        double dn = math::jacobi_dn(mu,u);
        double sn = sqrt(mu1) * sn_temp / dn;

        double sn_temp_plus = math::jacobi_sn(mu,u_plus);
        double dn_plus = math::jacobi_dn(mu,u_plus);
        double sn_plus = sqrt(mu1) * sn_temp_plus / dn_plus;


        theta = _pitchAmp * sn;

        theta_plus = _pitchAmp * sn_plus;


        thetadot = (theta_plus - theta) / _dt;

      }

      //
		return thetadot;
	}
	/////////////////////////////////////
  inline double gety(double time) const {
    /////////////////
    double u;

    double y;


    if (_m >= 0){
      u = time;

      // Jacobi function
      y = _plungeAmp * math::jacobi_sn(_m,u);
    }
    else {
      double mu = abs(_m) / (1 + abs(_m));
      double mu1 = 1 / (1 + abs(_m));
      u = time;


      double sn_temp = math::jacobi_sn(mu,u);
      double dn = math::jacobi_dn(mu,u);
      double sn = sqrt(mu1) * sn_temp / dn;

      y = _plungeAmp * sn;
    }

      //
    return y;
  }
  ///////////
    inline Motion* clone() const {

        return new PitchPlungeT_JF(
            _pitchAmp,

            _pitchPhase,
            _plungeAmp,

            _plungePhase,

            _dt,
            _m
        );
    };

private:
    double _pitchAmp;

    double _pitchPhase;
    double _plungeAmp;

    double _plungePhase;

    double _dt;
    double _m;
};

} // namespace ibpm

#endif /* _PITCHPLUNGET_JF_H_ */
