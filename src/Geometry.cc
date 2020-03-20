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
// $HeadURL$

#include "Geometry.h"
#include "Motion.h"
#include "RigidBody.h"
#include "Regularizer.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "utils.h"

using std::istringstream;
using std::getline;
using std::ifstream;

namespace ibpm {

Geometry::Geometry() {
    _numPoints = 0;
    _isStationary = true;
}

Geometry::Geometry(string filename) {
    _numPoints = 0;
    _isStationary = true;
    load( filename );
}

Geometry::~Geometry() {}

int Geometry::getNumPoints() const {
    return _numPoints;
}

Motion* Geometry::transferMotion() {
    /*
        This transfers the motion of the first RigidBody to an unsteady base flow.
        This is only implemented for a single RigidBody in motion.
        If there are more bodies, then we need to keep a copy of the motion 
          that we are about to delete, so that we can transform all of the 
          remaining bodies into a frame fixed to the first body.
        We go through the bodies, and take motion from the first one that is moving.  
    */
    vector<RigidBody>::iterator body;
    Motion* tempmotion = NULL;
    body = _bodies.begin();
    if( body->isStationary() == false ) {
        tempmotion = body->getMotion();
        body->clearMotion();  // this removes the motion and resets the isStationary flag to true.
    }
    // recompute whether or not remaining RigidBody objects are stationary
    _isStationary = true;
    for (body = _bodies.begin(); body != _bodies.end(); ++body) {
        if( body->isStationary() == false ) _isStationary = false;
    }
    return tempmotion->clone(); 
} 

void Geometry::transferCenter(double &x, double &y) {
    // This fills (x,y) with the center (xC,yC) from the first RigidBody. 
    vector<RigidBody>::iterator body;
    body = _bodies.begin();
    body->getCenter(x,y);
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
}

void Geometry::addBody(const RigidBody& body) {
    _bodies.push_back(body);
    _numPoints += body.getNumPoints();
    _isStationary = _isStationary && body.isStationary();
}

// Input format is as follows:
//    name My Geometry
//    body Flat Plate
//        line 0 0 1 0 0.1  # points on a line, spacing approximately 0.1
//        center 0.25 0     # center at quarter chord
//        motion fixed 0 0 0.3 # 0.3 radians angle of attack
//    end
//    body Large Circle
//        circle 2 3 5 0.1 # Points on a circle
//                         # default center is (2,3)
//    end
//    body Airfoil
//        raw naca0012.in  # Read in the raw data file
//    end
bool Geometry::load(istream& in) {
#ifdef DEBUG
    cout << "In Geometry::load(in), reading from input stream" << endl;
    cout << "  in.good() = " << in.good() << endl;
#endif
    string buf;
    string cmd;
    bool error_found = false;
    while ( getline( in, buf ) ) {
#ifdef DEBUG
        cerr << "You typed:" << buf << endl;
#endif
        istringstream one_line( buf );
        one_line >> cmd;
        MakeLowercase( cmd );
        if ( one_line.eof() || cmd[0] == '#' ) {
#ifdef DEBUG
            cerr << "[comment]" << endl;
#endif
        }
        else if ( cmd == "body" ) {
            string name;
            getline( one_line, name );
            EatWhitespace( name );
            RigidBody body;
            body.setName( name );
            if ( body.load( in ) ) {
                addBody( body );
#ifdef DEBUG
                cerr << "Load a new body, named [" << name << "]" << endl;
#endif
            }
            else error_found = true;
        }
        else if ( cmd == "end" ) {
#ifdef DEBUG
            cerr << "End" << endl;
#endif
            break;
        }
        else {
            cerr << "WARNING: could not parse the following line:" << endl;
            cerr << buf << endl;
            error_found = true;
        }
    }
    if ( error_found ) return false;
    else return true;
}

bool Geometry::load(string filename) {
#ifdef DEBUG
    cout << "in Geometry::load, opening file " << filename << endl;
#endif
    ifstream in( filename.c_str() );
    if ( in.good() ) {
        return load( in );
    }
    else {
        cerr << "Error: could not open " << filename << " for input." << endl;
        return false;
    }
}

} // namespace ibpm
