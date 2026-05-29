#ifndef _PITCHPLUNGEFS_H_
#define _PITCHPLUNGEFS_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

/*!
    \file PitchPlunge.h
    \class PitchPlunge

    \brief Subclass of Motion, for a pitching and plunging body

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 29 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class PitchPlungeFS : public Motion {
public:

    /// \brief Define a Motion corresponding to sinusoidal pitching and
    /// plunging, centered about the origin:
    ///    y(t)      = A_y      sin( 2*pi * f_y t + phi_y)
    ///    \theta(t) = A_\theta sin( 2*pi * f_theta t + phi_theta)
	PitchPlungeFS(
		//
		double pitchAmplitude0,
		double pitchWaveNum0,

		double pitchAmplitude1,
		double pitchWaveNum1,

		double pitchAmplitude2,
		double pitchWaveNum2,

		double pitchAmplitude3,
		double pitchWaveNum3,

		double pitchAmplitude4,
		double pitchWaveNum4,

		double pitchAmplitude5,
		double pitchWaveNum5,

		double pitchAmplitude6,
		double pitchWaveNum6,


		double Tt,

		double pitchPhase,
		double plungeWaveNum,
		double plungeAmplitude,
		double plungePhase

	) :
		_pitchAmp0(pitchAmplitude0),
		_pitchWN0(pitchWaveNum0),

		_pitchAmp1(pitchAmplitude1),
		_pitchWN1(pitchWaveNum1),

		_pitchAmp2(pitchAmplitude2),
		_pitchWN2(pitchWaveNum2),

		_pitchAmp3(pitchAmplitude3),
		_pitchWN3(pitchWaveNum3),

		_pitchAmp4(pitchAmplitude4),
		_pitchWN4(pitchWaveNum4),

		_pitchAmp5(pitchAmplitude5),
		_pitchWN5(pitchWaveNum5),

		_pitchAmp6(pitchAmplitude6),
		_pitchWN6(pitchWaveNum6),

		_Tp(Tt),
		_pitchPhase(pitchPhase),
		_plungeAmp(plungeAmplitude),
		_plungeWN(plungeWaveNum),
		_plungePhase(plungePhase)
		 {
		double twopi = 8. * atan(1.);
			_Tp = 1.0;
			_pitchWN0 /= (twopi);
			_pitchWN1 /= (twopi);
			_pitchWN2 /= (twopi);
			_pitchWN3 /= (twopi);
			_pitchWN4 /= (twopi);
			_pitchWN5 /= (twopi);
			_pitchWN6 /= (twopi);

	}


    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, y(t), theta(t), 0, ydot(t), thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
        // double y    = _plungeAmp * sin( _plungeFreq * time + _plungePhase );
        // double ydot = _plungeAmp * _plungeFreq * cos( _plungeFreq * time + _plungePhase );
        // double theta    = _pitchAmp * sin( _pitchFreq * time + _pitchPhase );
        // double thetadot = _pitchAmp * _pitchFreq * cos( _pitchFreq * time + _pitchPhase );

        //double y = 0;
        //double ydot = 0;
        //double theta = 0;
        //double thetadot = 0;
/*
		double y = _plungeAmp * sin(1 / _plungeWN * twopi * time / _Tp + _plungePhase);
         double ydot = _plungeAmp / _plungeWN * twopi / _Tp * cos( 1 / _plungeWN * twopi * time / _Tp + _plungePhase );*/

		double y = _plungeAmp * sin(1. / 1.  * time / _Tp + _plungePhase) * 0.;
		 double ydot = _plungeAmp / 1.  / _Tp * cos(1. / 1. * time / _Tp + _plungePhase) * 0.;

		 double theta0 = _pitchAmp0 * sin(1. / _pitchWN0 * time / _Tp + _pitchPhase);
         double thetadot0 = _pitchAmp0 / _pitchWN0 / _Tp * cos( 1. / _pitchWN0 * time / _Tp + _pitchPhase );

		 double theta1 = _pitchAmp1 * sin(1. / _pitchWN1 * time / _Tp + _pitchPhase);
		 double thetadot1 = _pitchAmp1 / _pitchWN1 / _Tp * cos(1. / _pitchWN1 * time / _Tp + _pitchPhase);

		 double theta2 = _pitchAmp2 * sin(1. / _pitchWN2 * time / _Tp + _pitchPhase);
		 double thetadot2 = _pitchAmp2 / _pitchWN2  / _Tp * cos(1. / _pitchWN2 * time / _Tp + _pitchPhase);

		 double theta3 = _pitchAmp3 * sin(1. / _pitchWN3 * time / _Tp + _pitchPhase);
		 double thetadot3 = _pitchAmp3 / _pitchWN3 / _Tp * cos(1. / _pitchWN3 * time / _Tp + _pitchPhase);

		 double theta4 = _pitchAmp4 * sin(1. / _pitchWN4 * time / _Tp + _pitchPhase);
		 double thetadot4 = _pitchAmp4 / _pitchWN4 / _Tp * cos(1. / _pitchWN4 * time / _Tp + _pitchPhase);

		 double theta5 = _pitchAmp5 * sin(1. / _pitchWN5 * time / _Tp + _pitchPhase);
		 double thetadot5 = _pitchAmp5 / _pitchWN5 / _Tp * cos(1. / _pitchWN5 * time / _Tp + _pitchPhase);

		 double theta6 = _pitchAmp6 * sin(1. / _pitchWN6 * time / _Tp + _pitchPhase);
		 double thetadot6 = _pitchAmp6 / _pitchWN6 / _Tp * cos(1. / _pitchWN6 * time / _Tp + _pitchPhase);

		 double theta = theta0 + theta1 + theta2 + theta3 + theta4 + theta5 + theta6;
		 double thetadot = thetadot0 + thetadot1 + thetadot2 + thetadot3 + thetadot4 + thetadot5 + thetadot6;

		 //double theta = theta0 + theta1 + theta2;
		 //double thetadot = thetadot0 + thetadot1 + thetadot2;

    cout << "time step..." << time << endl;
   //  cout << "Tp..." << _Tp << endl;
    cout << "theta..." << theta << endl;
		cout << "theta_dot..." << thetadot << endl;
	//cout << " time Fun... " << time << endl;
		 //cout << "()..." << 1 / _pitchWN0 * time / _Tp + _pitchPhase << endl;

		 //cout << "time..." << time << endl;
		 //cout << "Tp..." << _Tp << endl;
		 //cout << "pitchPahse..." << _pitchPhase << endl;
		 //cout << "amp0..." << _pitchAmp0<< endl;
		 cout << "WN1..." << _pitchWN1 << endl;
		 cout << "pitchAmp1..." << _pitchAmp1 << endl;

		 //cout << "y..." << y<< endl;
		 //cout << "ydot..." << ydot << endl;

        return TangentSE2( 0, y, theta, 0, ydot, thetadot );
    }

    inline Motion* clone() const {
		double twopi = 8. * atan(1.);
        return new PitchPlungeFS(
			_pitchAmp0,
			_pitchWN0 * (twopi),

			_pitchAmp1,
			_pitchWN1 * (twopi),

			_pitchAmp2,
			_pitchWN2 * (twopi),

			_pitchAmp3,
			_pitchWN3 * (twopi),

			_pitchAmp4,
			_pitchWN4 * (twopi),

			_pitchAmp5,
			_pitchWN5 * (twopi),

			_pitchAmp6,
			_pitchWN6 * (twopi),

			_Tp,
			_pitchPhase,
			_plungeAmp,
			_plungeWN,
			_plungePhase

        );
    };

private:
    double _pitchAmp0;
    double _pitchWN0;

	double _pitchAmp1;
	double _pitchWN1;

	double _pitchAmp2;
	double _pitchWN2;

	double _pitchAmp3;
	double _pitchWN3;

	double _pitchAmp4;
	double _pitchWN4;

	double _pitchAmp5;
	double _pitchWN5;

	double _pitchAmp6;
	double _pitchWN6;


	double _Tp;
    double _pitchPhase;
    double _plungeAmp;
    double _plungeWN;
    double _plungePhase;

};

} // namespace ibpm

#endif /* _PITCHPLUNGE_H_ */
