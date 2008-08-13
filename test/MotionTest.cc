#include "Motion.h"
#include "FixedPosition.h"
#include "TangentSE2.h"
#include <gtest/gtest.h>

TEST( Motion, FixedPosition ) {
    double x = 1;
    double y = 2;
    double theta = 3;

    Motion* motion;
    motion = new FixedPosition( x, y, theta );

    // Check stationary flag is true
    EXPECT_EQ( motion->isStationary(), true );

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
    
}