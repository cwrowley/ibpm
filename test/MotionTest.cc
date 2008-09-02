#include "FixedPosition.h"
#include "PitchPlunge.h"
#include "TangentSE2.h"
#include <gtest/gtest.h>

using namespace ibpm;

TEST( Motion, FixedPosition ) {
    double x = 1;
    double y = 2;
    double theta = 3;

    Motion* motion;
    motion = new FixedPosition( x, y, theta );

    // Check stationary flag is true
    EXPECT_EQ( true, motion->isStationary() );

    // Check position is correct
    TangentSE2 g = motion->getTransformation(0.);
    double a, b, omega;
    g.getPosition( a, b, omega );
    EXPECT_DOUBLE_EQ( a, x );
    EXPECT_DOUBLE_EQ( b, y );
    EXPECT_DOUBLE_EQ( omega, theta );    

    // Check velocity is zero
    double adot, bdot, omegadot;
    g.getVelocity( adot, bdot, omegadot );
    EXPECT_DOUBLE_EQ( adot, 0 );
    EXPECT_DOUBLE_EQ( bdot, 0 );
    EXPECT_DOUBLE_EQ( omegadot, 0 );
 
    delete motion;
}

TEST( Motion, PitchPlunge ) {
    double pitchAmp = 2;
    double plungeAmp = 3;
    double pitchFreq = 4;
    double plungeFreq = 5;
    const double twopi = 8. * atan(1.);
    
    Motion* motion;
    motion = new PitchPlunge( pitchAmp, pitchFreq, plungeAmp, plungeFreq );
    
    // Check stationary flag is false
    EXPECT_EQ( false, motion->isStationary() );
    
    // Check position is correct
    TangentSE2 g = motion->getTransformation(0.);
    double x, y, theta;
    g.getPosition( x, y, theta );
    EXPECT_DOUBLE_EQ( 0, x );
    EXPECT_DOUBLE_EQ( 0, y );
    EXPECT_DOUBLE_EQ( 0, theta );

    // Check velocity is correct
    double xdot, ydot, thetadot;
    // WORK HERE
    g.getVelocity( xdot, ydot, thetadot );
    EXPECT_DOUBLE_EQ( 0, xdot );
    EXPECT_DOUBLE_EQ( plungeAmp * twopi * plungeFreq, ydot );
    EXPECT_DOUBLE_EQ( pitchAmp * twopi * pitchFreq, thetadot );
    // TODO: add phase shift to pitch/plunge
    // TODO: check different points in time
}
