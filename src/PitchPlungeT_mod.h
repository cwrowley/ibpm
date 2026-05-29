#ifndef _PITCHPLUNGET_MOD_H_
#define _PITCHPLUNGET_MOD_H_

#include <array>
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

class PitchPlungeT_mod : public Motion {
public:

    /// \brief Define a Motion corresponding to sinusoidal pitching and
    /// plunging, centered about the origin:
    ///    y(t)      = A_y      sin( 2*pi * f_y t + phi_y)
    ///    \theta(t) = A_\theta sin( 2*pi * f_theta t + phi_theta)
    PitchPlungeT_mod(
        double dt,
        double *motiontime,
        double *theta_in,
        double *y_in,
        int size
        ) :
        _dt(dt),
        _motiontime(motiontime),
        _theta_in(theta_in),
        _y_in(y_in),
        _size(size) {

    }

    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, y(t), theta(t), 0, ydot(t), thetadot(t))
    inline TangentSE2 getTransformation(double time) const {
      ///////////////
        // if (time < 1){
        //   time = time + _size;
        // }
        //////////////
        double ind_temp = floor(time / _dt);
        int index = (int) ind_temp;
        if (index < 0){    ///   <0
          index = index + _size;
        }
        // double y = _y_in[index];

        double y;
        for (int i = 0; i<=_size; ++i){
          if (time >= (double)i * _dt && time < ((double)i + 1) * _dt){
            if (i >= _size){
              y = _y_in[0] + (_y_in[0+1]-_y_in[0])/_dt*(time-(double)index * _dt);

            }

            else{
              y = _y_in[index] + (_y_in[index+1]-_y_in[index])/_dt*(time-(double)index * _dt);
            }

          }
        }
        double ydot_i;
        if (index < 1){
          // ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);///////////// for nonperiodic signal
          ydot_i = (_y_in[_size] - _y_in[_size-1]) / (_dt);
        }
        else if(index >= _size){
          ydot_i = (_y_in[index] - _y_in[index-1]) / (_dt);
        }
        else{
          ydot_i = (_y_in[index] - _y_in[index-1]) / (_dt);
        }
        double ydot_plus;
        if (index < 1){
          // ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);///////////// for nonperiodic signal
          ydot_plus = (_y_in[index+1] - _y_in[index]) / (_dt);
        }
        else if(index >= _size){
          ydot_plus = (_y_in[1] - _y_in[0]) / (_dt);
        }
        else{
          ydot_plus = (_y_in[index+1] - _y_in[index]) / (_dt);
        }
        //
        double ydot;
        for (int i = 0; i<=_size; ++i){
          if (time >= (double)i * _dt && time < ((double)i + 1) * _dt){
             ydot = ydot_i + (ydot_plus-ydot_i)/_dt*(time-(double)index * _dt);
          }
        }




        //////////////


                double theta;
                for (int i = 0; i<=_size; ++i){
                  if (time >= (double)i * _dt && time < ((double)i + 1) * _dt){
                    if (i >= _size){
                      theta = _theta_in[0] + (_theta_in[0+1]-_theta_in[0])/_dt*(time-(double)index * _dt);

                    }
                    else{
                      theta = _theta_in[index] + (_theta_in[index+1]-_theta_in[index])/_dt*(time-(double)index * _dt);
                    }

                  }
                }
                double thetadot_i;
                if (index < 1){
                  // ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);/////////////
                  thetadot_i = (_theta_in[_size] - _theta_in[_size-1]) / (_dt);
                }
                else if(index >= _size){
                  thetadot_i = (_theta_in[index] - _theta_in[index-1]) / (_dt);
                }
                else{
                  thetadot_i = (_theta_in[index] - _theta_in[index-1]) / (_dt);
                }
                double thetadot_plus;
                if (index < 1){
                  // ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);/////////////
                  thetadot_plus = (_theta_in[index+1] - _theta_in[index]) / (_dt);
                }
                else if(index >= _size){
                  thetadot_plus = (_theta_in[1] - _theta_in[0]) / (_dt);
                }
                else{
                  thetadot_plus = (_theta_in[index+1] - _theta_in[index]) / (_dt);
                }
                //
                double thetadot;
                for (int i = 0; i<=_size; ++i){
                  if (time >= (double)i * _dt && time < ((double)i + 1) * _dt){
                     thetadot = thetadot_i + (thetadot_plus-thetadot_i)/_dt*(time-(double)index * _dt);
                  }
                }



        // cout << " theta_minus ... " << _theta_in[index-1] << endl;
        // cout << " theta ... " << _theta_in[index] << endl;
        // cout << " theta_plus ... " << _theta_in[index+1] << endl;
        // cout << " thetadot ... " << thetadot << endl;
		//cout << "ydot_true..." << ydot << endl;
        return TangentSE2( 0, y, theta, 0, ydot, thetadot );
		//                 x           x_dot   x_dot = 1 ?

    }
	///////////////////////////////////
	inline double getydot(double time) const {
    ///////////////
      // if (time < 1){
      //   time = time + _size;
      // }
      //////////////
    //////////////
    double ind_temp = floor(time / _dt);
    int index = (int) ind_temp;
    if (index < 0){
      index = index + _size;
    }
    double ydot_i;
    if (index < 1){
      // ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);/////////////
      ydot_i = (_y_in[_size] - _y_in[_size-1]) / (_dt);
    }
    else if(index >= _size){
      ydot_i = (_y_in[index] - _y_in[index-1]) / (_dt);
    }
    else{
      ydot_i = (_y_in[index] - _y_in[index-1]) / (_dt);
    }
    double ydot_plus;
    if (index < 1){
      // ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);/////////////
      ydot_plus = (_y_in[index+1] - _y_in[index]) / (_dt);
    }
    else if(index >= _size){
      ydot_plus = (_y_in[1] - _y_in[0]) / (_dt);
    }
    else{
      ydot_plus = (_y_in[index+1] - _y_in[index]) / (_dt);
    }
    //
    double ydot;
    for (int i = 0; i<=_size; ++i){
      if (time >= (double)i * _dt && time < ((double)i + 1) * _dt){
         ydot = ydot_i + (ydot_plus-ydot_i)/_dt*(time-(double)index * _dt);
      }
    }

    cout << " index ... " << index << endl;
    cout << " size ... " << _size << endl;
    cout << " y_minus ... " << _y_in[index-1] << endl;
    cout << " y ... " << _y_in[index] << endl;
    cout << " y_plus ... " << _y_in[index+1] << endl;
    cout << " time ... " << time << endl;
    cout << " ydot ... " << ydot << endl;
		return ydot;
	}

  inline double gety(double time) const {
    ///////////////
      // if (time < 1){
      //   time = time + _size;
      // }
      //////////////
    //////////////
    double ind_temp = floor(time / _dt);
    int index = (int) ind_temp;
    if (index < 0){
      index = index + _size;
    }
    double y;
    for (int i = 0; i<=_size; ++i){
      if (time >= (double)i * _dt && time < ((double)i + 1) * _dt){
        if (i >= _size){
          y = _y_in[0] + (_y_in[0+1]-_y_in[0])/_dt*(time-(double)index * _dt);
        }
        else{
          y = _y_in[index] + (_y_in[index+1]-_y_in[index])/_dt*(time-(double)index * _dt);
        }

      }
    }
    return y;
  }

	inline double gettheta(double time) const {
    ///////////////
      // if (time < 1){
      //   time = time + _size;
      // }
      //////////////
		// double theta = _pitchAmp * sin(time / _pitchT + _pitchPhase);
    //////////////
    double ind_temp = floor(time / _dt);
    int index = (int) ind_temp;
    if (index < 0){
      index = index + _size;
    }
    double theta;
    for (int i = 0; i<=_size; ++i){
      if (time >= (double)i * _dt && time < ((double)i + 1) * _dt){
        if (i >= _size){
          theta = _theta_in[0] + (_theta_in[0+1]-_theta_in[0])/_dt*(time-(double)index * _dt);
        }
        else{
          theta = _theta_in[index] + (_theta_in[index+1]-_theta_in[index])/_dt*(time-(double)index * _dt);
        }

      }
    }
		return theta;
	}

	inline double getthetadot(double time) const {
    ///////////////
      // if (time < 1){
      //   time = time + _size;
      // }
      //////////////
		// double thetadot = _pitchAmp / _pitchT * cos(time / _pitchT + _pitchPhase);
    //////////////
    double ind_temp = floor(time / _dt);
    int index = (int) ind_temp;
    if (index < 0){
      index = index + _size;
    }
    double thetadot_i;
    if (index < 1){
      // ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);/////////////
      thetadot_i = (_theta_in[_size] - _theta_in[_size-1]) / (_dt);
    }
    else if(index >= _size){
      thetadot_i = (_theta_in[index] - _theta_in[index-1]) / (_dt);
    }
    else{
      thetadot_i = (_theta_in[index] - _theta_in[index-1]) / (_dt);
    }
    double thetadot_plus;
    if (index < 1){
      // ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);/////////////
      thetadot_plus = (_theta_in[index+1] - _theta_in[index]) / (_dt);
    }
    else if(index >= _size){
      thetadot_plus = (_theta_in[1] - _theta_in[0]) / (_dt);
    }
    else{
      thetadot_plus = (_theta_in[index+1] - _theta_in[index]) / (_dt);
    }
    //
    double thetadot;
    for (int i = 0; i<=_size; ++i){
      if (time >= (double)i * _dt && time < ((double)i + 1) * _dt){
         thetadot = thetadot_i + (thetadot_plus-thetadot_i)/_dt*(time-(double)index * _dt);
      }
    }
		return thetadot;
	}

	/////////////////////////////////////
    inline Motion* clone() const {
        double twopi = 8. * atan(1.);
        return new PitchPlungeT_mod(
            _dt,
            _motiontime,
            _theta_in,
            _y_in,
            _size
        );
    };

private:
    double _dt;
    double *_motiontime;
    double *_theta_in;
    double *_y_in;
    int _size;

};

} // namespace ibpm

#endif /* _PITCHPLUNGET_H_ */
