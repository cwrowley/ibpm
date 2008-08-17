// Geometry.cc
//
// Description:
//
// Author(s):
// Clancy Rowley
//
// Date: 3 Jul 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL:// $Header$



#include "Geometry.h"
#include "RigidBody.h"

Geometry::Geometry() {
    _numPoints = 0;
}

Geometry::~Geometry() {}

int Geometry::getNumPoints() const {
    return _numPoints;
}

BoundaryVector Geometry::getPoints() const {
    BoundaryVector coords(_numPoints);
    vector<RigidBody>::const_iterator body;

    int ind = 0;
    for (body = _bodies.begin(); body != _bodies.end(); ++body) {
        BoundaryVector bodyCoords = body->getPoints();
        for (int bodyInd=0; bodyInd < bodyCoords.getNumPoints(); ++bodyInd) {
            coords(X,ind) = bodyCoords(X,bodyInd);
            coords(Y,ind) = bodyCoords(Y,bodyInd);
            ++ind;
        }
    }
    return coords;
}

BoundaryVector Geometry::getVelocities() const {
    BoundaryVector velocities(_numPoints);
    velocities = 0;
    return velocities;
}

bool Geometry::isStationary() const { 
    // TODO: call the isStationary methods of the associated RigidBodies
    // or save an instance variable with the value of this flag
    return true;
}

void Geometry::moveBodies(double time) const {
    // TODO: Update the positions of the associated RigidBodies and
    // recompute Regularizer operations
}

void Geometry::addBody(const RigidBody& body) {
    _bodies.push_back(body);
    _numPoints += body.getNumPoints();
}

void Geometry::load(const istream& in) {
    
}
