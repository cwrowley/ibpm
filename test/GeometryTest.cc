#include "RigidBody.h"
#include "BoundaryVector.h"
#include "Geometry.h"
#include <gtest/gtest.h>

namespace {

class GeometryTest : public testing::Test {
protected:
    GeometryTest() {
        // Define a RigidBody with one point
        RigidBody body1;
        body1.addPoint(0.,0.);
        
        // Define a RigidBody wth two points
        RigidBody body2;
        body2.addPoint(0., 1.);
        body2.addPoint(0., -1.);
        
        _geom.addBody(body1);
        _geom.addBody(body2);
    }

    // data
    Geometry _geom;
};

TEST_F(GeometryTest, NoBodies) {
    Geometry emptyGeom;
    
    EXPECT_EQ( emptyGeom.getNumPoints(), 0 );

    BoundaryVector points = emptyGeom.getPoints();
    EXPECT_EQ( points.getNumPoints(), 0 );
    
    BoundaryVector velocities = emptyGeom.getVelocities();
    EXPECT_EQ( velocities.getNumPoints(), 0 );    
}

TEST_F(GeometryTest, OneBody) {
    const double x = 3.;
    const double y = 5.;
    RigidBody body;
    body.addPoint(x,y);
    
    Geometry geom;
    EXPECT_EQ( geom.getNumPoints(), 0 );
    geom.addBody(body);    
    EXPECT_EQ( geom.getNumPoints(), 1 );
    
    BoundaryVector pts = geom.getPoints();
    EXPECT_EQ( pts.getNumPoints(), 1 );
    EXPECT_DOUBLE_EQ( pts(X,0), x );
    EXPECT_DOUBLE_EQ( pts(Y,0), y );
}

TEST_F(GeometryTest, TwoBodies) {
    // Coordinates of points in bodies to test
    const double x[3] = {3, 5, 7};
    const double y[3] = {10, 20, 30};

    // Define a RigidBody with one point
    RigidBody body1;
    body1.addPoint( x[0], y[0] );
    Geometry geom;
    geom.addBody(body1);
    EXPECT_EQ( geom.getNumPoints(), 1 );
    
    // Define a RigidBody wth two points
    RigidBody body2;
    body2.addPoint( x[1], y[1] );
    body2.addPoint( x[2], y[2] );
    geom.addBody(body2);
    EXPECT_EQ( geom.getNumPoints(), 3 );

    BoundaryVector points = geom.getPoints();
    for (int i=0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ( points(X,i), x[i] );
        EXPECT_DOUBLE_EQ( points(Y,i), y[i] );
    }
}


}  // namespace

