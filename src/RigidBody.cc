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

#include "Grid.h"
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

/// For now numPoints is determined by circumference/gridspacing
void RigidBody::addCircle(double xc, double yc, double radius, int numPoints) {
    double x,y;
    const double pi = 4. * atan(1.);
    double dx = .01;
    numPoints = floor(2*pi*radius/dx);  
    for(int i=0;i<numPoints;i++) {
       x = xc + radius*cos(2*pi*i/numPoints);
       y = yc + radius*sin(2*pi*i/numPoints);
       addPoint(x,y);
    }   
}

void RigidBody::addLine(double x1, double y1, double x2, double y2, int numPoints) {
    double x,y;
    double dx = .01;
    double length = sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));
    numPoints = floor(length/dx);
    for(int i=0;i<=numPoints;i++) {
       x = x1 + (i/numPoints)*(x2-x1);
       y = y1 + (i/numPoints)*(y2-y1);
       addPoint(x,y);
    }
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
