#ifndef _MOTION_H_
#define _MOTION_H_

class TangentSE2;

/*!
    \file Motion.h
    \class Motion

    \brief Specify a position and orientation as a function of time
    
    Abstract base class: must be subclassed to be instantiated

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 12 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class Motion {
public:
    /// True if the body is moving (default is False)
    virtual bool isStationary() { return false; }

    /// Return a Euclidean transformation and its velocity
    /// (an element of TSE(2) at the specified time
    virtual TangentSE2 getTransformation(double time) = 0;
};

#endif /* _MOTION_H_ */
