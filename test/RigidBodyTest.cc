#include "RigidBody.h"
#include "Direction.h"
#include "BoundaryVector.h"
#include <gtest/gtest.h>

double tol = 1e-14;

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

TEST_F( RigidBodyTest, AddLine ) {
    // Add a line with 101 points from (1,1) to (3,4)
    body.addLine(1., 1., 3., 4., 101);
    BoundaryVector b = body.getPoints();

    // Check that numPoints is correct
    EXPECT_EQ( body.getNumPoints(), 101 );

    // Check that all spaces are equal
    for( int i = 0; i < 99; i++ ) {
       EXPECT_NEAR( b(X,i+1) - b(X,i) , b(X,i+2) - b(X,i+1) , tol );
       EXPECT_NEAR( b(Y,i+1) - b(Y,i) , b(Y,i+2) - b(Y,i+1) , tol );
    }
   
    // Check that start, end, and center points are correct
    EXPECT_DOUBLE_EQ( b(X,0), 1.);
    EXPECT_DOUBLE_EQ( b(Y,0), 1.);
    EXPECT_DOUBLE_EQ( b(X,100), 3.);
    EXPECT_DOUBLE_EQ( b(Y,100), 4.);
    EXPECT_DOUBLE_EQ( b(X,50), 2.);
    EXPECT_DOUBLE_EQ( b(Y,50), 2.5);
}

TEST_F( RigidBodyTest, AddCircle ) {
    double xp, yp;
    double d1, d2;
    // Add a circle with 360 points w/ center (2.5,2.5) and radius 2.5
    body.addCircle( 2.5, 2.5, 2.5, 360 );
    BoundaryVector b = body.getPoints();

    // Check that numPoints is correct
    EXPECT_EQ( body.getNumPoints(), 360 );
 
    // Check that radius is correct for every point
    for( int i = 0; i < 360; i++ ) {
       xp = b(X,i) - 2.5;
       yp = b(Y,i) - 2.5;
       EXPECT_NEAR( sqrt( xp*xp + yp*yp ) , 2.5 , tol ); 
    }  

    // Check that spacing is the same each time
    for( int i = 0; i < 358; i++ ) {
       d1 = sqrt( (b(X,i+1)-b(X,i))*(b(X,i+1)-b(X,i)) + (b(Y,i+1)-b(Y,i))*(b(Y,i+1)-b(Y,i)) );      
       d2 = sqrt( (b(X,i+2)-b(X,i+1))*(b(X,i+2)-b(X,i+1)) + (b(Y,i+2)-b(Y,i+1))*(b(Y,i+2)-b(Y,i+1)) );
       EXPECT_NEAR( d1, d2, tol );
    }

    // Check that 0, pi/2, pi, 3*pi/2 are at the correct places
    EXPECT_NEAR( b(X, 0), 5., tol );
    EXPECT_NEAR( b(Y, 0), 2.5, tol ); 
    EXPECT_NEAR( b(X, 90), 2.5, tol );
    EXPECT_NEAR( b(Y, 90), 5., tol );
    EXPECT_NEAR( b(X, 180), 0., tol );
    EXPECT_NEAR( b(Y, 180), 2.5, tol );
    EXPECT_NEAR( b(X, 270), 2.5, tol );
    EXPECT_NEAR( b(Y, 270), 0., tol );

}
