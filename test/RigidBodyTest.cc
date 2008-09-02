#include "RigidBody.h"
#include "Direction.h"
#include "BoundaryVector.h"
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>

using namespace std;
using namespace ibpm;

double tol = 1e-14;
double lowtol = 1e-5;

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
    body.addLine_n(1., 1., 3., 4., 101);
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
    double dx, dy;
    double xc = 2.; // coordinates of center of circle
    double yc = 3.; 
    double radius = 1.5;
    int numPoints = 8;

    // Add a circle
    body.addCircle_n( xc, yc, radius, numPoints );
    BoundaryVector b = body.getPoints();

    // Check that numPoints is correct
    EXPECT_EQ( body.getNumPoints(), numPoints );
 
    // Check that radius is correct for every point
    for( int i = 0; i < numPoints; i++ ) {
        dx = b(X,i) - xc;
        dy = b(Y,i) - yc;
        EXPECT_NEAR( dx*dx + dy*dy , radius * radius , tol ); 
    }

    // Compute distance^2 between first two points
    dx = b(X,1) - b(X,0);
    dy = b(Y,1) - b(Y,0);
    double distanceSquared = dx * dx + dy * dy;

    // Check that spacing is the same each time
    for( int i = 0; i < numPoints; i++ ) {
        int iNext = (i+1) % numPoints;
        dx = b(X,iNext) - b(X,i);
        dy = b(Y,iNext) - b(Y,i);
        EXPECT_NEAR( distanceSquared, dx*dx + dy*dy, tol );
    }
    
    // find indices of points at north, south, east, and west points on circle
    int east = 0;
    int north = numPoints / 4;
    int west = numPoints / 2;
    int south = 3 * numPoints / 4;
    
    // Check that north, south, east, and west are at the correct places
    EXPECT_NEAR( b(X, east),  xc + radius, tol );
    EXPECT_NEAR( b(Y, east),  yc,          tol ); 
    EXPECT_NEAR( b(X, north), xc,          tol );
    EXPECT_NEAR( b(Y, north), yc + radius, tol );
    EXPECT_NEAR( b(X, west),  xc - radius, tol );
    EXPECT_NEAR( b(Y, west),  yc,          tol );
    EXPECT_NEAR( b(X, south), xc,          tol );
    EXPECT_NEAR( b(Y, south), yc - radius, tol );

}

TEST_F( RigidBodyTest, IORaw1 ) {
    double x1 = 1.;
    double y1 = 2.;
    double x2 = 1./3;
    double y2 = 1./4;
    double x3 = 4. * atan(1.);
    double y3 = sqrt(2.);
    filebuf fb;
    fb.open ("./output/saveRaw.dat",ios::out);
    ostream out(&fb);

    // Add Three Points
    body.addPoint(x1,y1);
    body.addPoint(x2,y2);
    body.addPoint(x3,y3);
    body.saveRaw(out);   
    fb.close();
}

// TEST_F( RigidBodyTest, IORaw2 ) {
//     int numPoints = 3;
//     double x1 = 1.;
//     double y1 = 2.;
//     double x2 = 1./3;
//     double y2 = 1./4;
//     double x3 = 4. * atan(1.);
//     double y3 = sqrt(2);
//     filebuf fb;
//     fb.open("saveRaw.tmp",ios::in);
//     istream in(&fb);
// 
//     body.loadRaw(in);
//     fb.close();
//     BoundaryVector b = body.getPoints();
// 
//     // Check that numPoints is correct
//     EXPECT_EQ( body.getNumPoints() , numPoints );
//     
//     // Check that three points are correct
//     EXPECT_NEAR( b(X, 0), x1 , lowtol );
//     EXPECT_NEAR( b(Y, 0), y1 , lowtol );
//     EXPECT_NEAR( b(X, 1), x2 , lowtol ); 
//     EXPECT_NEAR( b(Y, 1), y2 , lowtol );
//     EXPECT_NEAR( b(X, 2), x3 , lowtol );
//     EXPECT_NEAR( b(Y, 2), y3 , lowtol );
// }
