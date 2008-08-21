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
    return _points.size();
};

void RigidBody::saveRaw(ostream& out) {
    int n = getNumPoints();
    out << n;
    for(int i=0; i<n; i++) {
       out << "\n" << setw(10) << _points[i].x << setw(10) << _points[i].y; 	
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
