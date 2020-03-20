#ifndef _MOTIONFILE_H_
#define _MOTIONFILE_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "utils.h"

using std::ifstream;


namespace ibpm {

// Defining a struct for storing (t,x,y,theta), which I call TxTSE2
struct TxTSE2 {
    TxTSE2(double t_in, double x_in, double y_in, double theta_in, double dx_in, double dy_in, double dtheta_in) {
        t = t_in;
        x = x_in;
        y = y_in;
        theta = theta_in;
        dx = dx_in;
        dy = dy_in;
        dtheta = dtheta_in;
    }
    double t;
    double x;
    double y;
    double theta;
    double dx;
    double dy;   
    double dtheta;
};

/*!
    \file MotionFile.h
    \class MotionFile

    \brief Subclass of Motion, for a motions defined in a file

    \author Steven Brunton
    \author $LastChangedBy: sbrunton $
    \date 26 Apr 2010
    \date $LastChangedDate: 2008-09-03 02:41:03 -0400 (Wed, 03 Sep 2008) $
    \version $Revision: 132 $
*/

class MotionFile : public Motion {
public:
    
    /// \brief Define a Motion corresponding to data specified in a file
    MotionFile(
        string filename
        ) :
        _filename(filename) {
        ifstream in( filename.c_str() );
        int n;
        double t,x,y,theta,dx,dy,dtheta;
        double tlast = -1.e+5;
        in >> n;
        if( in.fail() ) {
            cerr << "ERROR:: MotionFile: file " << filename <<" formatted incorrectly! (err1)" << endl;
        }
        for(int i=0; i<n; i++) {
            in >> t >> x >> y >> theta >> dx >> dy >> dtheta;
            if( in.fail() ) {
                cerr << "ERROR:: MotionFile: file " << filename << " formatted incorrectly! (err2)" << endl;
            }
            if( t<tlast ) {
                cerr << "ERROR:: MotionFile: time must increase monotonically!" << endl;
            }
            tlast = t;
            addTxTSE2(t,x,y,theta,dx,dy,dtheta);
        }   
    }

    /// Adds an element of TxTSE2 to the vector
    void addTxTSE2(double t, double x, double y, double theta, double dx, double dy, double dtheta) {
        TxTSE2 p(t,x,y,theta,dx,dy,dtheta);
        // Add the element to the motion file data list
        _data.push_back(p);
    }
    
    /// Returns transformation for the piecewise linear data in _filename:
    inline TangentSE2 getTransformation(double time) const {
        int index = -1;
        // find correct index corresponding to current time
        int n=_data.size();
        for(int i=0; i<n-1; i++) {
            if( (time>=_data[i].t) && (time<_data[i+1].t) ) index = i;
        } 
        double xo = 0.;
        double yo = 0.;
        double thetao = 0.;
        double dxo = 0.;
        double dyo = 0.;
        double dthetao = 0.;
        double alpha=0.;
        double tdiff=0.;
        if( index != -1 ) {
            tdiff = _data[index+1].t-_data[index].t;
            alpha = (time-_data[index].t)/tdiff;
            xo = alpha*_data[index+1].x+(1-alpha)*_data[index].x;
            yo = alpha*_data[index+1].y+(1-alpha)*_data[index].y;
            thetao = alpha*_data[index+1].theta+(1-alpha)*_data[index].theta;
            dxo = alpha*_data[index+1].dx+(1-alpha)*_data[index].dx;
            dyo = alpha*_data[index+1].dy+(1-alpha)*_data[index].dy;
            dthetao = alpha*_data[index+1].dtheta+(1-alpha)*_data[index].dtheta;
        }
        return TangentSE2( xo, yo, thetao, dxo, dyo, dthetao );
    }

    inline Motion* clone() const {
        return new MotionFile(
            _filename
        );
    };

private:
    string _filename;
    vector<TxTSE2> _data;
};
} // namespace ibpm

#endif /* _MOTIONFILE_H_ */
