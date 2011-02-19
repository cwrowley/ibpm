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


using namespace std;
using namespace ibpm;

namespace ibpm {

BaseFlow::BaseFlow() {
    _xCenter = 0;
	_yCenter = 0;
	_isStationary = true;
	_motion = NULL;
    _time = 0.;
	_mag = 0.;
	_alpha = 0.;
}

BaseFlow::BaseFlow(const Grid& grid ) {
    _xCenter = 0;
    _yCenter = 0;
	_isStationary = true;
	_motion = NULL;
    _time = 0.;
    resize( grid );
	_mag = 0.;
	_alpha = 0.;
	_q = Flux::UniformFlow( grid, _mag, _alpha);
}

BaseFlow::BaseFlow(const Grid& grid, double mag, double alpha) {
    _xCenter = 0;
	_yCenter = 0;
	_isStationary = true;
	_motion = NULL;
	_time = 0.;
	_mag = mag;
	_alpha = alpha;
	resize( grid );
	_q = Flux::UniformFlow( grid, _mag, _alpha);
}

BaseFlow::BaseFlow(const Grid& grid, double mag, double alpha, const Motion& motion) {
    _xCenter = 0;
	_yCenter = 0;
	_motion = motion.clone();
	_isStationary = motion.isStationary();
	_motion = NULL;
	_time = 0.;
	_mag = mag;
	_alpha = alpha;
	resize( grid );
	_q = Flux::UniformFlow( grid, _mag, _alpha);
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

void BaseFlow::moveFlow(double time) {
    if ( _motion == NULL ) return;

        double x,y,theta,xdot,ydot,thetadot;
        double xdotnew, ydotnew;
	TangentSE2 g = _motion->getTransformation(time);
        g.getPosition(x,y,theta);
        g.getVelocity(xdot,ydot,thetadot);
 
        // The true xdot and ydot are determined by g\in TSE2 and the constant base flow (alpha,mag)
        xdotnew = _mag*cos(theta+_alpha) - xdot;
        ydotnew = _mag*sin(theta+_alpha) - ydot;
        TangentSE2 gnew(x,y,theta,xdotnew,ydotnew,thetadot);
        cout << "xdotnew = " << xdotnew << endl;
        cout << "ydotnew = " << ydotnew << endl;
        cout << "thetadot = " << thetadot << endl;
        cout << "x = " << x << " , y = " << y << " , th = " << theta << " alpha = " << _alpha << " , mag = " << _mag << endl;
        // Update the baseFlow based on this new motion 
        Grid grid = _q.getGrid();
        _q.TSE2Flow( grid, gnew, _xCenter, _yCenter);	
}

BaseFlow::~BaseFlow() {}
    
} // namespace ibpm
