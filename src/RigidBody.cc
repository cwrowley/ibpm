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
#include "FixedPosition.h"
#include "FixedVelocity.h"
#include "PitchPlunge.h"
#include "SigmoidalStep.h"
#include "LagStep1.h"
#include "LagStep2.h"
#include "EldredgeManeuver.h"
#include "EldredgeCombined2.h"
#include "Eldredge1.h"
#include "Eldredge2.h"
#include "MotionFile.h"
#include "MotionFilePeriodic.h"
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "utils.h"

using std::setw;
using std::ifstream;
using std::istringstream;

namespace ibpm {

RigidBody::RigidBody() {
    _xCenter = 0;
    _yCenter = 0;
    _isStationary = true;
    _motion = NULL;
}

RigidBody::RigidBody(const RigidBody& body) :
   _xCenter( body._xCenter ),
   _yCenter( body._yCenter ),
   _isStationary( body._isStationary ),
   _refPoints( body._refPoints ),
   _currentPoints( body._currentPoints ),
   _currentVelocities( body._currentVelocities ) {
   if ( body._motion == NULL ) {
       _motion = NULL;
   }
   else {
       _motion = body._motion->clone();
   }
}

RigidBody::~RigidBody() {
    if ( _motion != NULL ) {
        cout << "         deleting _motion" << endl;
        delete _motion;
    }
}

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
	// To round a value x, take floor( x + 0.5 )
    int numPoints = (int) floor( twopi / dTheta + 1 + 0.5 );
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
	// To round a value x, take floor( x + 0.5 )
    int numPoints = (int) floor( length / dx + 1 + 0.5 );
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

void RigidBody::addLine_aoa(
    double l,
    double xC,
    double yC,
    double aoa,
    int numPoints
    ) {
        double x0, y0, x0r, y0r, x1r, y1r, x1, y1;
	double pi =  3.141592653589793238462643383279502884197169399375;
	double alpha = aoa / 180 * pi;
        double cosa = cos(alpha);
        double sina = sin(alpha);
        double delta = l / (numPoints - 1);
        for(int i=0; i < numPoints; i++) {
            x0 = i * delta;
            y0 = 0.;    // (x0,y0) point before rotation
            x0r = x0-xC;
            y0r = y0-yC; // (x0r,y0r) referenced to center (xC,yC)
            x1r = x0r*cosa - y0r*sina;
            y1r = x0r*sina + y0r*cosa; // (x1r,y1r) after rotation, referenced to center (xC,yC)
            x1 = x1r+xC;
            y1 = y1r+yC;
            addPoint(x1,y1);
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
//     motion fixed x y theta
//     motion pitchplunge amp1 freq1 amp2 freq2
//     end
// Whitespace at the beginning of the line is ignored
// Returns false if invalid input was encountered
bool RigidBody::load(istream& in) {
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
        else if ( cmd == "line_aoa" ) {
            double l;
            double xC;
            double yC;
            double aoa;
            int numPoints;
            one_line >> l >> xC >> yC >> aoa >> numPoints;
            RB_CHECK_FOR_ERRORS;
            addLine_aoa( l, xC, yC, aoa, numPoints );
            setCenter( xC, yC );
#ifdef DEBUG
            cerr << "Add a line: length" << l << ", AoA" << aoa  
                ", n = " << numPoints << endl;
#endif
        }
        else if ( cmd == "motion" ) {
            string motionType;
            one_line >> motionType;
            MakeLowercase( motionType );
            if ( motionType == "fixed" ) {
                // FixedPosition
                double x;
                double y;
                double theta;
                one_line >> x >> y >> theta;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new FixedPosition( x, y, theta );
                setMotion( *m );
            }
            else if ( motionType == "fixedvel" ) {
                // FixedVelocity
                double xdot;
                double ydot;
                double thetadot;
                one_line >> xdot >> ydot >> thetadot;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new FixedVelocity( xdot, ydot, thetadot );
                setMotion( *m );
            }
            else if ( motionType == "pitchplunge" ) {
                // PitchPlunge
                double amp1;
                double freq1;
                double phase1;
                double amp2;
                double freq2;
                double phase2;
                one_line >> amp1 >> freq1 >> phase1 >> amp2 >> freq2 >> phase2;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new PitchPlunge( amp1, freq1, phase1, amp2, freq2, phase2 );
                setMotion( *m );
            }
            else if ( motionType == "sigmoidalstep" ) {
                // SigmoidalStep
                double AMP;
                double DUR;
                double startTime;
                one_line >> AMP >> DUR >> startTime;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new SigmoidalStep( AMP, DUR, startTime);
                setMotion( *m ); 
            }
            else if ( motionType == "lagstep1" ) {
                // First-order lag filtered square pulse
                double AMP;
                double PW;
                double TAU;
                double T0;
                one_line >> AMP >> PW >> TAU >> T0;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new LagStep1( AMP, PW, TAU, T0);
                setMotion( *m );
            }
            else if ( motionType == "lagstep2" ) {
                // Second-order lag filtered square pulse
                double AMP;
                double PW;
                double TAU;
                double T0;
                one_line >> AMP >> PW >> TAU >> T0;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new LagStep2( AMP, PW, TAU, T0);
                setMotion( *m );
            }
            else if ( motionType == "eldredge" ) {
                // Eldredge's canonical maneuver
                double AMP;
                double a;
                double t1;
                double t2;
                double t3;
                double t4;
                one_line >> AMP >> a >> t1 >> t2 >> t3 >> t4;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new EldredgeManeuver( AMP, a, t1, t2, t3, t4);
                setMotion( *m );
            }
            else if ( motionType == "eldredgecombined2" ) {
                // combining type eldredge for pitch and eldredge2 for plunge
                double AMPa;
                double a;
                double a1;
                double a2;
                double a3;
                double a4;
                double AMPb;
                double b;
                double b1;
                double b2;
                double b3;
                double b4;
                one_line >> AMPa >> a >> a1 >> a2 >> a3 >> a4 >> AMPb >> b >> b1 >> b2 >> b3 >> b4;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new EldredgeCombined2( AMPa, a, a1, a2, a3, a4, AMPb, b, b1, b2, b3, b4);
                setMotion( *m );
            }
            else if ( motionType == "eldredge1" ) {
                // Eldredge's canonical maneuver for plunging
                double AMP;
                double a;
                double t1;
                double t2;
                double t3;
                double t4;
                one_line >> AMP >> a >> t1 >> t2 >> t3 >> t4;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new Eldredge1( AMP, a, t1, t2, t3, t4);
                setMotion( *m );
            }
            else if ( motionType == "eldredge2" ) {
                // Eldredge's canonical maneuver for plunging
                double AMP;
                double a;
                double t1;
                double t2;
                double t3;
                double t4;
                one_line >> AMP >> a >> t1 >> t2 >> t3 >> t4;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new Eldredge2( AMP, a, t1, t2, t3, t4);
                setMotion( *m );
            }
            else if ( motionType == "motionfile" ) {
                // Motion defined in a file
                string filename;
                one_line >> filename;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new MotionFile( filename );
                setMotion( *m );
            }
            else if ( motionType == "motionfileperiodic" ) {
                // Periodic motion defined in a file
                string filename;
                double period;
                one_line >> filename >> period;
                RB_CHECK_FOR_ERRORS;
                Motion* m = new MotionFilePeriodic( filename, period );
                setMotion( *m );
            }
        }
        else if ( cmd == "name" ) {
            string name;
            getline( one_line, name );
            RB_CHECK_FOR_ERRORS;
            EatWhitespace( name );
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
    // Delete the old motion, if present
    if ( _motion != NULL ) {
        delete _motion;
    }
    // make a local copy of the new motion
    _motion = motion.clone();
    _isStationary = motion.isStationary();
}

Motion* RigidBody::getMotion() {
    return _motion->clone();
}

void RigidBody::clearMotion() {
    _motion = NULL;
    _isStationary = true;
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
        g.mapPosition(p->x-_xCenter, p->y-_yCenter, xnew, ynew);
        Point q(xnew+_xCenter,ynew+_yCenter);
        _currentPoints.push_back(q);

        // add the new velocity to the list of current velocities
        double u;
        double v;
        g.mapVelocity(p->x-_xCenter, p->y-_yCenter, u, v);
        Point vel(u,v);
        _currentVelocities.push_back(vel);
    }
}

} // namespace ibpm
