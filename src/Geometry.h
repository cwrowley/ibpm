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
	// ...
};

/*! \brief Test routine for now
\param[in] i Integer argument
\return argument incremented by 1
*/
int plus_one(int i);

#endif /* _GEOMETRY_H_ */
