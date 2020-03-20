// BaseFlow.cc
//
// Description:
// Implementation of the BaseFlow class
//
// Author(s):
// Steve Brunton
//
// Date: 9 Nov 2009
//
// $Revision: 288 $
// $LastChangedDate: 2009-11-09 14:58:56 -0500 (Mon, 09 Nov 2009) $
// $LastChangedBy: sbrunton $
// $HeadURL: svn+ssh://sbrunton@rainier.princeton.edu/ibpm/branches/ibpm-unsteadyBaseflow/src/BaseFlow.cc $

#include "BaseFlow.h"
#include "Grid.h"
#include "Direction.h"
#include "BoundaryVector.h"
#include "TangentSE2.h"
#include "Motion.h"
#include "FixedPosition.h"
#include "PitchPlunge.h"
#include "Flux.h"
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "utils.h"


using namespace ibpm;

namespace ibpm {

BaseFlow::BaseFlow() {
    _xCenter = 0;
    _yCenter = 0;
    _isStationary = true;
    _motion = NULL;
    _time = 0.;
    _mag = 0.;
    _magBF = 0.;
    _alpha = 0.;
    _alphaBF = 0.;
    _gamma = 0.;
}

BaseFlow::BaseFlow(const Grid& grid ) {
    _xCenter = 0;
    _yCenter = 0;
    _isStationary = true;
    _motion = NULL;
    _time = 0.;
    resize( grid );
    _mag = 0.;
    _magBF = 0.;
    _alpha = 0.;
    _alphaBF = 0.;
    _gamma = 0.;
    _q = Flux::UniformFlow( grid, _magBF, _alphaBF);
}

BaseFlow::BaseFlow(const Grid& grid, double mag, double alpha) {
    _xCenter = 0;
    _yCenter = 0;
    _isStationary = true;
    _motion = NULL;
    _time = 0.;
    _mag = mag;
    _magBF = mag;
    _alpha = alpha;
    _alphaBF = alpha;
    _gamma = -1.*_alphaBF;
    resize( grid );
    _q = Flux::UniformFlow( grid, _magBF, _alphaBF);
}

BaseFlow::BaseFlow(const Grid& grid, double mag, double alpha, const Motion& motion) {
    _xCenter = 0;
    _yCenter = 0;
    _motion = motion.clone();
    _isStationary = motion.isStationary();
    _motion = NULL;
    _time = 0.;
    _mag = mag;
    _magBF = mag;
    _alpha = alpha;
    _alphaBF = alpha;
    _gamma = -1.*_alphaBF;
    resize( grid );
    _q = Flux::UniformFlow( grid, _magBF, _alphaBF);
}

void BaseFlow::resize( const Grid& grid ) {
    _q.resize( grid );
}

bool BaseFlow::isStationary() const {
    return _isStationary;
}

void BaseFlow::setMotion(const Motion& motion) {
    // Delete the old motion, if present
    if ( _motion != NULL ) {
        delete _motion;
    }
    // make a local copy of the new motion
    _motion = motion.clone();
    _isStationary = motion.isStationary();
}

void BaseFlow::setAlphaMag(double time) {
    double x,y,theta,xdot,ydot,thetadot;   /// Motion of rigid body
    double xdotBF, ydotBF;    /// Velocity components of base flow (Uinf,alphaBF)
    double xdotT, ydotT;      /// Total velocity (sum of xdot and xdotBF)
    TangentSE2 g = _motion->getTransformation(time);
    g.getPosition(x,y,theta);
    g.getVelocity(xdot,ydot,thetadot);
    xdotBF = _magBF*cos(_alphaBF);
    ydotBF = _magBF*sin(_alphaBF);
    xdotT = xdot - xdotBF;
    ydotT = ydot - ydotBF;
    _gamma = atan2(ydotT,-1.*xdotT);
    _alpha = -1.*theta - _gamma;
    _mag = sqrt( xdotT*xdotT + ydotT*ydotT );
}

void BaseFlow::moveFlow(double time) {
    if ( _motion == NULL ) return;
    double x,y,theta,xdot,ydot,thetadot;   /// Motion of rigid body
    TangentSE2 g = _motion->getTransformation(time);
    g.getPosition(x,y,theta);
    g.getVelocity(xdot,ydot,thetadot);
    /* The flow is decomposed into a base flow of magnitude _mag at
        an angle _alpha = -theta-_gamma, and a purely rotational component
        $-\dot{\theta}$ centered at the body center of rotation. 
        The formulae are:
          _mag = (_magBF*cos(_alphaBF) - \dot{x}, -_magBF*sin(_alphaBF) + \dot{y})
          _alpha = -theta - gamma
    */
    setAlphaMag(time);

    xdot = _mag * cos(_alpha);  // using the xdot,ydot of the baseflow
    ydot = _mag * sin(_alpha);
    TangentSE2 gnew(x,y,theta,xdot,ydot,-1.*thetadot);
    // Update the baseFlow based on this new motion 
    _q.setFlow( gnew, _xCenter, _yCenter);	
}

BaseFlow::~BaseFlow() {}
    
} // namespace ibpm
