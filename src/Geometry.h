#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "BoundaryVector.h"
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
    int getNumPoints() const {
        return 1;
    }

    BoundaryVector getPoints() const {
        BoundaryVector coords(1);
        coords = 0;
        return coords;
    }
    
    BoundaryVector getVelocities() const {
        BoundaryVector velocities(1);
        velocities = 0;
        return velocities;
    }

    inline bool isStationary() const { 
        // TODO: call the isStationary methods of the associated RigidBodies
        // or save an instance variable with the value of this flag
        return true;
    }

    inline void moveBodies(double time) const {
        // TODO: Update the positions of the associated RigidBodies and
        // recompute Regularizer operations
    }
    
private:
    
};

/*! \brief Test routine for now
\param[in] i Integer argument
\return argument incremented by 1
*/
int plus_one(int i);

#endif /* _GEOMETRY_H_ */
