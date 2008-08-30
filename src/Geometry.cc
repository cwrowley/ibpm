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
#include "Regularizer.h"

namespace ibpm {

Geometry::Geometry() {
    _numPoints = 0;
    _isStationary = true;
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
    vector<RigidBody>::const_iterator body;

    int ind = 0;
    for (body = _bodies.begin(); body != _bodies.end(); ++body) {
        BoundaryVector bodyVel = body->getVelocities();
        for (int bodyInd=0; bodyInd < bodyVel.getNumPoints(); ++bodyInd) {
            velocities(X,ind) = bodyVel(X,bodyInd);
            velocities(Y,ind) = bodyVel(Y,bodyInd);
            ++ind;
        }
    }
    return velocities;
}

bool Geometry::isStationary() const { 
    return _isStationary;
}

void Geometry::moveBodies(double time) const {
    // Update the positions and velocities of the associated RigidBodies
    vector<RigidBody>::const_iterator body;
    for (body = _bodies.begin(); body != _bodies.end(); ++body) {
        body->moveBody(time);
    }

    // recompute Regularizer operations
    if (_regularizer != NULL) {
        _regularizer->update();
    }
}

void Geometry::addBody(const RigidBody& body) {
    _bodies.push_back(body);
    _numPoints += body.getNumPoints();
    _isStationary = _isStationary && body.isStationary();
}

void Geometry::load(const istream& in) {
    // TODO: implement this
}

void Geometry::setRegularizer(Regularizer& reg) const {
    _regularizer = &reg;
}

} // namespace ibpm
