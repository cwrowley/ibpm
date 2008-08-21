#ifndef _RIGIDBODY_H_
#define _RIGIDBODY_H_

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

// using std::istream;
// using std::ostream;
// using std::strin

class BoundaryVector;
class Motion;
struct Point;

/*!
\file RigidBody.h
\class RigidBody

\brief Specify coordinates and center of a rigid body, and its motion.

Detailed description here...

\author Clancy Rowley
\author $LastChangedBy$
\date 12 Aug 2008
\date $LastChangedDate$
\version $Revision$
*/

class RigidBody {
public:
    /// Default constructor, initializes the center to (0,0).
    RigidBody();

    /// Destructor
    ~RigidBody();
    
    /// Add the specified point to the list of points on the body's boundary
    void addPoint(double x, double y);    
    
    /// Add a circle with center (xc, yc) and the given radius
    /// with the specified number of points.
    void addCircle(
        double xc,
        double yc,
        double radius,
        int numPoints
    );

    /// Add a line connecting (x1,y1) and (x2,y2)
    /// with the specified number of points
    void addLine(
        double x1,
        double y1,
        double x2,
        double y2,
        int numPoints
    );
    
    /// Load a list of commands from the specified input stream
    void load(istream& in);
    
    /// Load a list of points, in ASCII format, with one point per line.
    /// Assumes the center is (0,0)
    void loadRaw(istream& in);
    
    /// Save a list of points to the specified output stream
    void saveRaw(ostream& out);
    
    /// Return the number of points on the body's boundary
    int getNumPoints() const;
    
    /// Return the list of coordinates for each point on the body
    BoundaryVector getPoints() const;
    
    /// Return true if the body is not moving in time
    bool isStationary();

    /// Set the evolution of the current body (which may be stationary or not)
    void setMotion(const Motion& motion);
    
    /// Set the center of the body, about which rotations are defined
    inline void setCenter(double x, double y) {
        _xCenter = x;
        _yCenter = y;
    }

    /// Get the center of the body, about which rotations are defined
    inline void getCenter(double& x, double& y) const {
        x = _xCenter;
        y = _yCenter;
    }
    
    /// Set the name of the body
    void setName(string name);
    
    /// Return the name of the body
    string getName();

private:
    string _name;
    double _xCenter;  ///< x-coordinate of center
    double _yCenter;  ///< y-coordinate of center
    vector<Point> _points;
};

// Define a small class for keeping track of points in 2d
struct Point {
    Point(double x_in, double y_in) {
        x = x_in;
        y = y_in;
    }
    double x;
    double y;
};

#endif /* _RIGIDBODY_H_ */
