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
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <ctype.h>

using namespace std;

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
    double dx
    ) {
    double dTheta = dx / radius;
    double twopi = 8. * atan(1.);
    int numPoints = twopi / dTheta + 1;
    addCircle_n( xc, yc, radius, numPoints );
}

void RigidBody::addCircle_n(
    double xc,
    double yc,
    double radius,
    int numPoints
    ) {
    double twopi = 8. * atan(1.);
    double dTheta = twopi / numPoints;
    for(int i=0; i < numPoints; i++) {
        double x = xc + radius * cos( i * dTheta );
        double y = yc + radius * sin( i * dTheta );
        addPoint( x, y );
    }   
}

void RigidBody::addLine(
    double x1,
    double y1,
    double x2,
    double y2,
    double dx
    ) {
    double length = sqrt( (x2-x1) * (x2-x1) + (y2-y1) * (y2-y1) );
    int numPoints = length / dx + 1;
    addLine_n( x1, y1, x2, y2, numPoints );
}

void RigidBody::addLine_n(
    double x1,
    double y1,
    double x2,
    double y2,
    int numPoints
    ) {
    double deltaX = (x2 - x1) / (numPoints - 1);
    double deltaY = (y2 - y1) / (numPoints - 1);
    for(int i=0; i < numPoints; i++) {
        double x = x1 + i * deltaX;
        double y = y1 + i * deltaY;
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


bool RigidBody::loadRaw(string filename) {
    ifstream in( filename.c_str() );
    int n;
    double x,y;
    in >> n;
    if ( in.fail() ) {
        return false;
    }
    for(int i=0; i<n; i++) { 
       in >> x >> y;
       if ( in.fail() ) {
           return false;
       }
       addPoint(x,y); 
    }
    return true;
}

static void parse_error(const string& buf) {
    cerr << "WARNING: could not parse the following line:" << endl;
    cerr << buf << endl;
}

static bool check_bad_input(
    const istringstream& linestream,
    const string& line
    ) {
    if ( linestream.fail() ) {
        parse_error( line );
        return true;
    }
    else {
        return false;
    }
}

// 
static void eat_whitespace( string& s ) {
    while ( isblank( s[0] ) ) {
        s.erase(0,1);
    }
}

static void make_lowercase( string& s ) {
    for (unsigned int i=0; i<s.length(); ++i) {
        s[i] = tolower(s[i]);
    }
}

#define RB_CHECK_FOR_ERRORS                     \
    if ( check_bad_input( one_line, buf ) ) {   \
        error_found = true;                     \
        continue;                               \
    }

// Load a list of commands from the specified input stream
// Input format is as follows:
//     name Name of this object
//     center x y  # location of the center of the object
//     point x y   # add a point at this location
//     point x y
//     point x y
//     line x1 y1 x2 y2 dx
//     line_n x1 y1 x2 y2 npts
//     circle xc yc radius dx
//     circle_n xc yc radius npts
//     raw naca0012.dat
//     end
// Whitespace at the beginning of the line is ignored
// Returns false if invalid input was encountered
bool RigidBody::load(istream& in) {
    string buf;
    string cmd;
    bool error_found = false;
    while ( getline( cin, buf ) ) {
#ifdef DEBUG
        cerr << "You typed:" << buf << endl;
#endif
        istringstream one_line( buf );
        one_line >> cmd;
        make_lowercase( cmd );
        if ( cmd[0] == '#' ) {
#ifdef DEBUG
            cerr << "[comment]" << endl;
#endif
        }
        else if ( cmd == "center" ) {
            double x;
            double y;
            one_line >> x >> y;
            RB_CHECK_FOR_ERRORS;
            setCenter( x, y );
#ifdef DEBUG
            cerr << "Specify the center: (" << x << ", " << y << ")" << endl;
#endif
        }
        else if ( cmd == "circle" ) {
            double xc;
            double yc;
            double radius;
            double dx;
            one_line >> xc >> yc >> radius >> dx;
            RB_CHECK_FOR_ERRORS;
            addCircle( xc, yc, radius, dx );
#ifdef DEBUG
            cerr << "Add a circle, center (" << xc << ", " << yc << "), "
                << "radius " << radius << ", dx = " << dx << endl;
#endif
        }
        else if ( cmd == "circle_n" ) {
            double xc;
            double yc;
            double radius;
            int numPoints;
            one_line >> xc >> yc >> radius >> numPoints;
            RB_CHECK_FOR_ERRORS;
            addCircle_n( xc, yc, radius, numPoints );
#ifdef DEBUG
            cerr << "Add a circle_n, center (" << xc << ", " << yc << "), "
                << "radius " << radius << ", n = " << numPoints << endl;
#endif
        }
        else if ( cmd == "end" ) {
#ifdef DEBUG
            cerr << "End" << endl;
#endif
            break;
        }
        else if ( cmd == "line" ) {
            double x0;
            double y0;
            double x1;
            double y1;
            double dx;
            one_line >> x0 >> y0 >> x1 >> y1 >> dx;
            RB_CHECK_FOR_ERRORS;
            addLine( x0, y0, x1, y1, dx );
#ifdef DEBUG
            cerr << "Add a line: (" << x0 << ", " << y0 << ") - ("
                << x1 << ", " << y1 << "), dx = " << dx << endl;
#endif
        }
        else if ( cmd == "line_n" ) {
            double x0;
            double y0;
            double x1;
            double y1;
            int numPoints;
            one_line >> x0 >> y0 >> x1 >> y1 >> numPoints;
            RB_CHECK_FOR_ERRORS;
            addLine_n( x0, y0, x1, y1, numPoints );
#ifdef DEBUG
            cerr << "Add a line: (" << x0 << ", " << y0 << ") - ("
                << x1 << ", " << y1 << "), n = " << numPoints << endl;
#endif
        }
        else if ( cmd == "name" ) {
            string name;
            getline( one_line, name );
            RB_CHECK_FOR_ERRORS;
            eat_whitespace( name );
            setName( name );
#ifdef DEBUG
            cerr << "Specify the name: " << name << endl;
#endif
        }
        else if ( cmd == "point" ) {
            double x;
            double y;
            one_line >> x >> y;
            RB_CHECK_FOR_ERRORS;
            addPoint( x, y );
#ifdef DEBUG
            cerr << "Add a point: (" << x << ", " << y << ")" << endl;            
#endif
        }
        else if ( cmd == "raw" ) {
            string filename;
            one_line >> filename;
            RB_CHECK_FOR_ERRORS;
            loadRaw( filename );
#ifdef DEBUG
            cerr << "Read a raw file: " << filename << endl;            
#endif
        }
        else {
            parse_error( buf );
        }
    }
    if ( error_found ) return false;
    else return true;
}

#undef RB_CHECK_FOR_ERRORS

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

void RigidBody::setName(string name) {
    _name = name;
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
