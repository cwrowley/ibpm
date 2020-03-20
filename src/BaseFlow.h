#ifndef _BASEFLOW_H_
#define _BASEFLOW_H_

#include "Grid.h"
#include "Flux.h"
#include "Scalar.h"
#include "BoundaryVector.h"
#include <string>

using std::string;

namespace ibpm {

class BoundaryVector;
class Motion;

/*!
    \file BaseFlow.h
    \class BaseFlow

    \brief Structure for a BaseFlow (flux) that moves with time

    \author Steve Brunton
    \author $LastChangedBy: sbrunton $
    \date  7 Jul 2008
    \date $LastChangedDate: 2009-11-09 14:58:56 -0500 (Mon, 09 Nov 2009) $
    \version $Revision: 288 $
*/

class BaseFlow {
public:
    /// Default constructor: do not allocate memory
    BaseFlow();
    
    /// Initialize without any constant or moving base flow
    BaseFlow( const Grid& grid );
	
    /// Initialize without any motion
    BaseFlow( const Grid& grid, double mag, double alpha);
	
    /// Initialize without any constant base flow
    BaseFlow( const Grid& grid, const Motion& motion );
	
    /// Initialize with constant and moving component
    BaseFlow( const Grid& grid, double mag, double alpha, const Motion& motion );

    ~BaseFlow();
    
    /// \brief Allocate memory, with the specified Grid and number of
    /// boundary points
    void resize( const Grid& grid );

    /// Return true if the body is not moving in time
    bool isStationary() const;

    /// Set the evolution of the current BaseFlow (which may be stationary or not)
    void setMotion(const Motion& motion);

    /// Set the center of the domain, about which rotations are defined
    inline void setCenter(double x, double y) {
        _xCenter = x;
        _yCenter = y;
    }

    /// Get the center of the domain, about which rotations are defined
    inline void getCenter(double& x, double& y) const {
        x = _xCenter;
        y = _yCenter;
    }

    /// Get the angle of attack 
    inline double getAlpha() const {
        return _alpha;
    }    

    /// Get the magnitude of base flow
    inline double getMag() const {
        return _mag;
    }

    /// Determine the magnitude and angle of the base flow, including rigid body motion
    void setAlphaMag(double time);

    /// Update the BaseFlow, based on the Motion
    void moveFlow(double time);
 
    /// Set the value of the flux _q
    inline void setFlux(Flux& f) { _q = f;};	

    inline void setFlux(double f) { _q = f;};

    /// Return flux _q
    inline const Flux& getFlux() const {return _q;};

private:
    double _magBF;    /// Initial magnitude of base flow	
    double _alphaBF;  /// Initial angle of base flow
    double _mag;      /// True magnitude of flow, taking into account rigid body motion
    double _gamma;    /// Flight Path Angle
    double _alpha;    /// True angle of attack, taking into account rigid body motion
    Flux _q;
    double _time;
    double _xCenter;
    double _yCenter;
    bool _isStationary;
    Motion *_motion;
};


} // namespace ibpm

#endif /* _BASEFLOW_H_ */
