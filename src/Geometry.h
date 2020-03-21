#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "BoundaryVector.h"
#include "Motion.h"
#include "RigidBody.h"
#include <iostream>
#include <vector>

using std::string;
using std::vector;

namespace ibpm {

class RigidBody;
class Regularizer;

/*!
    \file Geometry.h
    \class Geometry
    
    \brief Create, load, and save geometries composed of RigidBody objects.
    
    \author Clancy Rowley
    \date  3 Jul 2008
    
    $Revision$
    $LastChangedDate$
    $LastChangedBy$
    $HeadURL$
*/

class Geometry {
public:
    /// Constructor
    Geometry();
    /// Constructor: load the geometry from the specified file
    Geometry(string filename);
    
    /// Destructor
    ~Geometry();

    /// \brief Append the given RigidBody to the list of bodies in the 
    /// current geometry.
    /// Makes a copy of it internally.
    void addBody(const RigidBody& body);

    /// \brief Return number of boundary points
    int getNumPoints() const;

    /// \brief Return number of bodies
    inline int getNumBodies() const {
        return _bodies.size();
    }

    /// \brief Return motion from one of the bodies, and clear that motion
    /// useful for unsteady baseflow, since we want to initialize baseFlow's motion
    /// from the motion in a geometry.  We need to get that motion, and then remove it from the 
    /// rigid body object it came from.  
    Motion* transferMotion(); 

    /// \brief Fill (x,y) with the center of rotation of the first RigidBody
    void transferCenter(double &x, double &y);

    /// \brief Return the boundary points in the geometry
    BoundaryVector getPoints() const;
    
    /// \brief Return the velocities of the boundary points
    BoundaryVector getVelocities() const;

    /// \brief Return true if the body is not moving; false otherwise
    bool isStationary() const;

    /// \brief Move the boundary points and update their velocities
    void moveBodies(double time) const;

    /// \brief Load a geometry from the specified input stream.
    /// Returns false if invalid input was encountered
    /// Input format is as follows:
    ///    name My Geometry
    ///    body Flat Plate
    ///        line 0 0 1 0 0.1  # points on a line, spacing approximately 0.1
    ///        center 0.25 0     # center at quarter chord
    ///        motion fixed 0 0 0.3 # 0.3 radians angle of attack
    ///    end
    ///    body Large Circle
    ///        circle 2 3 5 0.1 # Points on a circle
    ///                         # default center is (2,3)
    ///    end
    ///    body Airfoil
    ///        raw naca0012.in  # Read in the raw data file
    ///    end
    bool load(istream& in);
    
    /// \brief Load a geometry from the specified input file
    /// Returns true if successful
    bool load(string filename);
    
private:
    vector<RigidBody> _bodies;
    int _numPoints;
    bool _isStationary;
};

} // namespace

#endif /* _GEOMETRY_H_ */
