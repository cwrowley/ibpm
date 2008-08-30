#ifndef _PITCHPLUNGE_H_
#define _PITCHPLUNGE_H_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

/*!
    \file PitchPlunge.h
    \class PitchPlunge

    \brief Subclass of Motion, for a pitching and plunging body

    \author Clancy Rowley
    \author $LastChangedBy: $
    \date 29 Aug 2008
    \date $LastChangedDate: $
    \version $Revision: $
*/

class PitchPlunge : public Motion {
public:
    
    /// \brief Define a Motion corresponding to sinusoidal pitching and 
    /// plunging, centered about the origin:
    ///    y(t)      = A_y      sin( 2*pi * f_y t)
    ///    \theta(t) = A_\theta sin( 2*pi * f_theta t)
    PitchPlunge(
        double pitchAmplitude,
        double pitchFrequency,
        double plungeAmplitude,
        double plungeFrequency
        ) :
        _pitchAmp(pitchAmplitude),
        _pitchFreq(pitchFrequency),
        _plungeAmp(plungeAmplitude),
        _plungeFreq(plungeFrequency) {
        
        double twopi = 8. * atan(1.);
        _pitchFreq *= twopi;
        _plungeFreq *= twopi;
    }
    
    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, y(t), theta(t), 0, ydot(t), thetadot(t))
    TangentSE2 getTransformation(double time) const {
        double y    = _plungeAmp * sin( _plungeFreq * time );
        double ydot = _plungeAmp * _plungeFreq * cos( _plungeFreq * time );
        double theta    = _pitchAmp * sin( _pitchFreq * time );
        double thetadot = _pitchAmp * _pitchFreq * cos( _pitchFreq * time );
        return TangentSE2( 0, y, theta, 0, ydot, thetadot );
    }

private:
    double _pitchAmp;
    double _pitchFreq;
    double _plungeAmp;
    double _plungeFreq;
};

} // namespace ibpm

#endif /* _PITCHPLUNGE_H_ */
