#ifndef _MOTIONGENERATOR_H_
#define _MOTIONGENERATOR_H_

#include <Eigen>
#include <Dense>
#include <Core>


/*!
    \file MotionGenerator.h
    \class MotionGenerator

    \brief Class for creating the motion profile given the parameters
    of Fourier modes.

*/

// The column 0 is pitching and column 1 is heaving
// The motion takes the form of
// y = sum (An * sin(2*pi*t/T + Bn))
// theta = sum (Cn * sin(2*pi*t/T + Dn))
// for Fourier mode = 0
// y = sum (An * sin(2*pi*t/T) + Bn * cos(2*pi*t/T))
// theta = sum (Cn * sin(2*pi*t/T) + Dn * cos(2*pi*t/T))
// for Fourier mode  = 1


// The input to the system is An, Bn, Cn, Dn, T, Fourier mode, t_length, motion_type, dt, PI
// t_length: length of te time
// motion_type: pitching, heaving or combined
// dt time step
// PI

using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

Mat MotionGenerator( Vec& A, Vec& B, Vec& C, Vec& D, double Fourier_mode, Vec& tau, int t_length, int motion_type, double dt, double PI){


  int size2 = motion_type;
  int size3 = A.size();
  Mat theta_temp(t_length + 1, size3);
  Mat y_temp(t_length + 1, size3);
  Mat MotionParameter(t_length + 1, size2+1);
  Vec MotionTime(t_length + 1);


  for (int k = 0; k < size3; ++k){
    for (int i = 0; i <= t_length; ++i){
        // switch between two sypes of Fourier modes
        if (Fourier_mode < 0.5){
          // y = sum (An * sin(t + Bn))
          // theta = sum (Cn * sin(t + Dn))
          theta_temp(i, k) =  C(k)* sin(2 * PI * (2*k+1) * i * dt/(tau(0)) + D(k));
          y_temp(i, k) =  A(k)* sin(2 * PI * (2*k+1) * i * dt/(tau(1)) + B(k));
          //
        }
        else if (Fourier_mode >0.5){
          // y = sum (An * sin + Bn * cos)
          // theta = sum (Cn * sin + Dn * cos)
          theta_temp(i, k) = C(k) * sin(2 * PI * (2*k + 1) * i * dt / (tau(0))) + D(k) * cos(2 * PI * (2*k + 1) * i * dt / (tau(0)));
          y_temp(i, k) = A(k) * sin(2 * PI * (2*k + 1) * i * dt / (tau(0))) + B(k) * cos(2 * PI * (2*k + 1) * i * dt / (tau(0)));
        }

          MotionTime(i) = i * dt;
  ////
          cout << "theta ... " << theta_temp(i, k) << endl;
          cout << "h ... " << y_temp(i, k) << endl;
          cout << "time ... " << MotionTime(i) << endl;
          cout << "k ... " << k << endl;
          cout << "D(k) ... " << D(k) << endl;
          cout << "C(k) ... " << C(k) << endl;
          cout << "tau(k) ... " << tau(0) << endl;
          cout << "sin ... " << sin(2 * PI * (2*k+1) * i * dt/(tau(0)) + D(k)) << endl;
          cout << "motion ... " << C(k)* sin(2 * PI * (2*k+1) * i * dt/(tau(0)) + D(k)) << endl;

    }


  }
 
  for (int j = 0; j < size2; ++j){
    for(int i = 0; i <= t_length; ++i){
      MotionParameter(i, j) = 0;
    }
  }
  
  for (int i = 0; i <= t_length; ++i) {
    for (int k = 0; k < size3; ++k) {
      MotionParameter(i, 0) = theta_temp(i, k) + MotionParameter(i, 0);
      MotionParameter(i, 1) = y_temp(i, k) + MotionParameter(i, 1);
    }
    MotionParameter(i, 2) = MotionTime(i);
  }
return MotionParameter;

}




#endif /* _BC_H_ */
