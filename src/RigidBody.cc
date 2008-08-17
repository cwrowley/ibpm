// RigidBody.cc
//
// Description:
// Implementation of the RigidBody class
//
// Author(s):
// Clancy Rowley
//
// Date: 12 Aug 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL:// $Header$

#include "Direction.h"
#include "BoundaryVector.h"
#include "RigidBody.h"


RigidBody::RigidBody() {
    _xCenter = 0;
    _yCenter = 0;
}

RigidBody::~RigidBody() {}

void RigidBody::addPoint(double x, double y) {
    Point p(x,y);
    _points.push_back(p);
}

int RigidBody::getNumPoints() const {
    return _points.size();
};

/// Return the list of coordinates for each point on the body
BoundaryVector RigidBody::getPoints() const {
    int n = getNumPoints();
    BoundaryVector pointList(n);

    for (int i=0; i<n; ++i) {
        pointList(X,i) = _points[i].x;
        pointList(Y,i) = _points[i].y;
    }
    return pointList;
}
