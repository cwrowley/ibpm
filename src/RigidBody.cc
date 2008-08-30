// RigidBody.cc
//
// Description:
// Implementation of the RigidBody class
//
// The points on the body are stored with respect to a reference configuration
// in _refPoints
// 
// The current locations of the points on the body, defined by the associated 
// Motion, are contained in _currentPoints, and are updated whenever 
// moveBodies() is called.
//
// Author(s):
// Clancy Rowley
//
// Date: 12 Aug 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "Grid.h"
#include "Direction.h"
#include "BoundaryVector.h"
#include "RigidBody.h"
#include "TangentSE2.h"
#include "Motion.h"

namespace ibpm {

RigidBody::RigidBody() {
    _xCenter = 0;
    _yCenter = 0;
    _isStationary = true;
    _motion = NULL;
}

RigidBody::~RigidBody() {}

void RigidBody::addPoint(double x, double y) {
    Point p(x,y);
    // Add the point to both lists: reference locations and current locations
    _refPoints.push_back(p);
    _currentPoints.push_back(p);
    Point zero(0,0);
    _currentVelocities.push_back(zero);
}

void RigidBody::addCircle(
    double xc,
    double yc,
    double radius,
    int numPoints
    ) {
    double x,y;
    const double pi = 4. * atan(1.);
    double dTheta = 2. * pi / numPoints;
    for(int i=0; i<numPoints; i++) {
        x = xc + radius * cos( i * dTheta );
        y = yc + radius * sin( i * dTheta );
        addPoint(x,y);
    }   
}

void RigidBody::addLine(
    double x1,
    double y1,
    double x2,
    double y2,
    int numPoints
    ) {
    double x,y;
    double deltaX = (x2 - x1) / (numPoints - 1);
    double deltaY = (y2 - y1) / (numPoints - 1);
    for(int i=0; i<numPoints; i++) {
        x = x1 + i * deltaX;
        y = y1 + i * deltaY;
        addPoint(x,y);
    }
} 

int RigidBody::getNumPoints() const {
    return _refPoints.size();
};

void RigidBody::saveRaw(ostream& out) {
    int n = getNumPoints();
    out << n;
    for(int i=0; i<n; i++) {
       out << "\n" << setw(10) << _refPoints[i].x
           << setw(10) << _refPoints[i].y; 	
    }
}


void RigidBody::loadRaw(istream& in) {
    int n = getNumPoints();
    double x,y;
    in >> n;
    for(int i=0; i<n; i++) { 
       in >> x >> y;
       addPoint(x,y); 
    }
}

string RigidBody::getName() {
    return _name;
}
    
BoundaryVector RigidBody::toBoundaryVector(const vector<Point> list) {
    int n = list.size();
    BoundaryVector BVList(n);

    for (int i=0; i<n; ++i) {
        BVList(X,i) = list[i].x;
        BVList(Y,i) = list[i].y;
    }
    return BVList;
}

/// Return the list of coordinates for each point on the body
BoundaryVector RigidBody::getPoints() const {
    return toBoundaryVector( _currentPoints );
}

/// Return the list of velocities at each point on the body
BoundaryVector RigidBody::getVelocities() const {
    return toBoundaryVector( _currentVelocities );
}

bool RigidBody::isStationary() const {
    return _isStationary;
}

void RigidBody::setMotion(const Motion& motion) {
    _motion = &motion;
    _isStationary = motion.isStationary();
}

// Update the position of the body (_currentPoints) based on the motion
void RigidBody::moveBody(double time) const {
    if ( _motion == NULL ) return;
    TangentSE2 g = _motion->getTransformation(time);
    _currentPoints.clear();
    _currentVelocities.clear();

    // for each reference point
    vector<Point>::const_iterator p;
    for (p = _refPoints.begin(); p != _refPoints.end(); ++p) {
        // add the new position to the list of current positions
        double xnew;
        double ynew;
        g.mapPosition(p->x, p->y, xnew, ynew);
        Point q(xnew,ynew);
        _currentPoints.push_back(q);

        // add the new velocity to the list of current velocities
        double u;
        double v;
        g.mapVelocity(p->x, p->y, u, v);
        Point vel(u,v);
        _currentVelocities.push_back(vel);
    }
}

} // namespace ibpm
