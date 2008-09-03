#ifndef _FIXEDPOSITION_H_
#define _FIXEDPOSITION_H_

#include "Motion.h"
#include "TangentSE2.h"

namespace ibpm {

/*!
    \file FixedPosition.h
    \class FixedPosition

    \brief Subclass of Motion, for stationary body
    
    Specifies the location of the center of the body, and an angle of
    rotation about the center

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 12 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class FixedPosition : public Motion {
public:
    
    /// Define a Motion corresponding to a fixed position and rotation
    /// about the center
    FixedPosition(
        double x,
        double y,
        double theta) :
        _x(x),
        _y(y),
        _theta(theta) {    
    }
    
    inline bool isStationary() const { return true; }
    
    /// Returns a fixed transformation for all time: (x,y,theta,0,0,0)
    inline TangentSE2 getTransformation(double time) const {
        return TangentSE2( _x, _y, _theta, 0, 0, 0 );
    }
    
    inline Motion* clone() const {
        return new FixedPosition( _x, _y, _theta );
    };

private:
    double _x;
    double _y;
    double _theta;
};

} // namespace ibpm

#endif /* _FIXEDPOSITION_H_ */
