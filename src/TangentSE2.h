#ifndef _TANGENTSE2_H_
#define _TANGENTSE2_H_

#include <math.h>

namespace ibpm {

/*!
    \file TangentSE2.h
    \class TangentSE2

    \brief An abstraction of the group TSE(2) of transformations in 2D
    
    Keeps track of translations and rotations in 2D, and their velocities.
    \author Clancy Rowley
    \author $LastChangedBy$
    \date 12 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class TangentSE2 {
public:
    
    /// Constructor: specify the base point and corresponding velocity
    /// The base point is (x,y,theta), velocity (xdot, ydot, thetadot)
    TangentSE2(
        double x,
        double y,
        double theta,
        double xdot,
        double ydot,
        double thetadot) :
        _x(x),
        _y(y),
        _theta(theta),
        _xdot(xdot),
        _ydot(ydot),
        _thetadot(thetadot) {
    }
    
    inline void setPosition(
        double x,
        double y,
        double theta
        ) {
        _x = x;
        _y = y;
        _theta = theta;    
    }

    inline void getPosition(
        double& x,
        double& y,
        double& theta
        ) {
        x = _x;
        y = _y;
        theta = _theta;
    }

    inline void setVelocity(
        double xdot,
        double ydot,
        double thetadot
        ) {
        _xdot = xdot;
        _ydot = ydot;
        _thetadot = thetadot;    
    }

    inline void getVelocity(
        double& xdot,
        double& ydot,
        double& thetadot
        ) {
        xdot = _xdot;
        ydot = _ydot;
        thetadot = _thetadot;
    }

    /// Given the point (a,b), compute the mapped point (a_new, b_new)
    inline void mapPosition(
        double a,
        double b,
        double& a_new,
        double& b_new
        ) {
        double cost = cos(_theta);
        double sint = sin(_theta);    
        a_new = _x + a * cost - b * sint;
        b_new = _y + a * sint + b * cost;
    }

    /// Given the point (a,b) (with zero initial velocity), compute the
    /// mapped velocity (u_new, v_new)
    inline void mapVelocity(
        double a,
        double b,
        double& u_new,
        double& v_new
        ) {
        double cost = cos(_theta);
        double sint = sin(_theta);
        u_new = - a * sint * _thetadot - b * cost * _thetadot + _xdot;
        v_new =   a * cost * _thetadot - b * sint * _thetadot + _ydot;
    }

private:
    // Note: may be better to implement this in terms of a vector (x,y,theta)
    // so that we can use matrix multiplication for transformations
    // Here, just the simplest implementation
    double _x;
    double _y;
    double _theta;
    double _xdot;
    double _ydot;
    double _thetadot;
};

}

#endif /* _TANGENTSE2_H_ */
