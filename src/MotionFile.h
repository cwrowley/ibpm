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


namespace ibpm {

// Defining a struct for storing (t,x,y,theta), which I call TxSE2
struct TxSE2 {
    TxSE2(double t_in, double x_in, double y_in, double theta_in) {
        t = t_in;
        x = x_in;
        y = y_in;
        theta = theta_in;
    }
    double t;
    double x;
    double y;
    double theta;
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
        double t,x,y,theta;
        double tlast = -1.e+5;
        in >> n;
        if( in.fail() ) {
            cerr << "ERROR:: MotionFile: file " << filename <<" formatted incorrectly! (err1)" << endl;
        }
        for(int i=0; i<n; i++) {
            in >> t >> x >> y >> theta;
            if( in.fail() ) {
                cerr << "ERROR:: MotionFile: file " << filename << " formatted incorrectly! (err2)" << endl;
            }
            if( t<tlast ) {
                cerr << "ERROR:: MotionFile: time must increase monotonically!" << endl;
            }
            tlast = t;
            addTxSE2(t,x,y,theta);
        }   
    }

    /// Adds an element of TxSE2 to the vector
    void addTxSE2(double t, double x, double y, double theta) {
        TxSE2 p(t,x,y,theta);
        // Add the element to the motion file data list
        _data.push_back(p);
    }
    
    /// Returns transformation for the piecewise linear data in _filename:
    inline TangentSE2 getTransformation(double time) const {
        int index = -1;
        // find correct index corresponding to current time
        int n=_data.size();
        for(int i=0; i<n-1; i++) {
            if( (time>_data[i].t) && (time<_data[i+1].t) ) index = i;
        }     
        double xout = 0.;
        double yout = 0.;
        double thetaout = 0.;
        double xdotout = 0.;
        double ydotout = 0.;
        double thetadotout = 0.;
        double alpha;
        double tdiff;
        if( index != -1 ) {
            tdiff = _data[index+1].t-_data[index].t;
            alpha = (time-_data[index].t)/tdiff;
            xout = alpha*_data[index+1].x+(1-alpha)*_data[index].x;
            yout = alpha*_data[index+1].y+(1-alpha)*_data[index].y;
            thetaout = alpha*_data[index+1].theta+(1-alpha)*_data[index].theta;
            xdotout = (_data[index+1].x-_data[index].x)/tdiff;
            ydotout = (_data[index+1].y-_data[index].y)/tdiff;
            thetadotout = (_data[index+1].theta-_data[index].theta)/tdiff;
        }
        return TangentSE2( xout, yout, thetaout, xdotout, ydotout, thetadotout );
    }

    inline Motion* clone() const {
        return new MotionFile(
            _filename
        );
    };

private:
    string _filename;
    vector<TxSE2> _data;
    vector<double[4]> _data2;
};
} // namespace ibpm

#endif /* _MOTIONFILE_H_ */
