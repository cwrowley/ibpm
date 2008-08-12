#include "RigidBody.h"
#include "Direction.h"
#include "BoundaryVector.h"
#include <gtest/gtest.h>

class RigidBodyTest : public testing::Test {
protected:
    RigidBody body;
};

TEST_F( RigidBodyTest, AddPoint ) {
    EXPECT_EQ( body.getNumPoints(), 0 );

    // Add one point
    body.addPoint(1.,2.);
    EXPECT_EQ( body.getNumPoints(), 1 );

    // Get the point in the form of a BoundaryVector and make sure it's equal
    BoundaryVector b = body.getPoints();
    EXPECT_EQ( b.getNumPoints(), 1 );
    EXPECT_DOUBLE_EQ( b(X,0), 1. );
    EXPECT_DOUBLE_EQ( b(Y,0), 2. );
    
    // Add another point
    body.addPoint(3.,4.);
    EXPECT_EQ( body.getNumPoints(), 2 );
    BoundaryVector b2 = body.getPoints();
    EXPECT_EQ( b2.getNumPoints(), 2 );
    EXPECT_DOUBLE_EQ( b2(X,0), 1. );
    EXPECT_DOUBLE_EQ( b2(Y,0), 2. );
    EXPECT_DOUBLE_EQ( b2(X,1), 3. );
    EXPECT_DOUBLE_EQ( b2(Y,1), 4. );
}

TEST_F( RigidBodyTest, SetCenter ) {
    // Add a point
    body.addPoint(1.,2.);
    
    // Check that the center is at the default (0,0)
    double xc;
    double yc;
    body.getCenter( xc, yc );
    EXPECT_DOUBLE_EQ( xc, 0. );
    EXPECT_DOUBLE_EQ( yc, 0. );

    // Change the center
    body.setCenter( 5., 6. );
    body.getCenter( xc, yc );
    EXPECT_DOUBLE_EQ( xc, 5. );
    EXPECT_DOUBLE_EQ( yc, 6. );
}