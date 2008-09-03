#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "BoundaryVector.h"
#include "RigidBody.h"
#include <iostream>
#include <vector>
using namespace std;

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
    
    /// Destructor
    ~Geometry();

    /// \brief Append the given RigidBody to the list of bodies in the 
    /// current geometry.
    /// Makes a copy of it internally.
    void addBody(const RigidBody& body);

    int getNumPoints() const;

    BoundaryVector getPoints() const;
    
    BoundaryVector getVelocities() const;

    bool isStationary() const;

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

    /// \brief Associate a Regularizer object.
    /// WARNING: Does not make a copy internally
    void setRegularizer(Regularizer& reg) const;
    
private:
    vector<RigidBody> _bodies;
    int _numPoints;
    bool _isStationary;
    mutable Regularizer* _regularizer;
};

} // namespace

#endif /* _GEOMETRY_H_ */
