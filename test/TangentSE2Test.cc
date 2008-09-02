#include "TangentSE2.h"
#include <gtest/gtest.h>
#include <math.h>

using namespace ibpm;

namespace {
    
const double pi = 4 * atan(1);

TEST( TangentSE2Test, Identity ) {
    TangentSE2 a( 0, 0, 0, 0, 0, 0 );
    double x = 1;
    double y = 2;

    double x1, y1;
    a.mapPosition( x, y, x1, y1 );
    EXPECT_DOUBLE_EQ( x, x1 );
    EXPECT_DOUBLE_EQ( y, y1 );
    
    double u, v;
    a.mapVelocity( x, y, u, v );
    EXPECT_DOUBLE_EQ( u, 0 );
    EXPECT_DOUBLE_EQ( v, 0 );
}

TEST( TangentSE2Test, Stationary ) {
    TangentSE2 a( 1, 2, 0, 0, 0, 0 );
    double x = 3;
    double y = 4;

    double x1, y1;
    a.mapPosition( x, y, x1, y1 );
    EXPECT_DOUBLE_EQ( x1, x + 1 );
    EXPECT_DOUBLE_EQ( y1, y + 2 );

    double u, v;
    a.mapVelocity( x, y, u, v );
    EXPECT_DOUBLE_EQ( u, 0 );
    EXPECT_DOUBLE_EQ( v, 0 );
}

TEST( TangentSE2Test, Rotate90 ) {
    TangentSE2 a( 0, 0, pi/2., 0, 0, 0 );
    double x = 2;
    double y = 3;
    
    double x1, y1;
    a.mapPosition( x, y, x1, y1 );
    EXPECT_DOUBLE_EQ( x1, -y );
    EXPECT_DOUBLE_EQ( y1, x );
}

TEST( TangentSE2Test, Moving ) {
    double u = 3;
    double v = 7;
    double omega = 11;
    TangentSE2 a( 0, 0, 0, u, v, omega );
    double x = 5;
    double y = 0;
    
    double xdot, ydot;
    a.mapVelocity( x, y, xdot, ydot );
    EXPECT_DOUBLE_EQ( xdot, u );
    EXPECT_DOUBLE_EQ( ydot, v + x * omega );    
}

TEST( TangentSE2Test, MovingWithTranslation ) {
    double u = 3;
    double v = 7;
    double omega = 11;
    TangentSE2 a( 1, 2, pi/2, u, v, omega );
    double x = 5;
    double y = 0;
    
    double xdot, ydot;
    a.mapVelocity( x, y, xdot, ydot );
    EXPECT_DOUBLE_EQ( xdot, u + -x * omega );
    EXPECT_DOUBLE_EQ( ydot, v + 0 );
}

} // namespace
