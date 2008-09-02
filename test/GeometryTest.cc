#include "RigidBody.h"
#include "BoundaryVector.h"
#include "Geometry.h"
#include "PitchPlunge.h"
#include <gtest/gtest.h>
#include <math.h>

using namespace ibpm;

namespace {

const double PI = 4. * atan(1.);

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
    int _nPoints;
    Geometry _geom;
};

#define EXPECT_ALL_EQ( a, b, nPoints )      \
    for (int i=0; i < nPoints; ++i) {       \
        EXPECT_DOUBLE_EQ( (a), (b) );       \
    }

TEST_F(GeometryTest, NoBodies) {
    Geometry emptyGeom;
    
    EXPECT_EQ( 0, emptyGeom.getNumPoints() );

    BoundaryVector points = emptyGeom.getPoints();
    EXPECT_EQ( 0, points.getNumPoints() );
    
    BoundaryVector velocities = emptyGeom.getVelocities();
    EXPECT_EQ( 0, velocities.getNumPoints() );    
}

TEST_F(GeometryTest, OneBody) {
    const double x = 3.;
    const double y = 5.;
    RigidBody body;
    body.addPoint(x,y);
    
    Geometry geom;
    EXPECT_EQ( 0, geom.getNumPoints() );
    geom.addBody(body);    
    EXPECT_EQ( 1, geom.getNumPoints() );
    
    BoundaryVector pts = geom.getPoints();
    EXPECT_EQ( pts.getNumPoints(), 1 );
    EXPECT_DOUBLE_EQ( pts(X,0), x );
    EXPECT_DOUBLE_EQ( pts(Y,0), y );
}

TEST_F(GeometryTest, TwoBodies) {
    // Coordinates of points in bodies to test
    const double x[3] = {3, 5, 7};
    const double y[3] = {10, 20, 30};

    Geometry geom;

    // Define a RigidBody with one point
    RigidBody* body;
    body = new RigidBody();
    body->addPoint( x[0], y[0] );
    geom.addBody( *body );
    delete body;
    EXPECT_EQ( 1, geom.getNumPoints() );

    // Define a RigidBody with two points
    body = new RigidBody();
    body->addPoint( x[1], y[1] );
    body->addPoint( x[2], y[2] );
    geom.addBody( *body );
    delete body;
    EXPECT_EQ( 3, geom.getNumPoints() );

    BoundaryVector points = geom.getPoints();
    for (int i=0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ( points(X,i), x[i] );
        EXPECT_DOUBLE_EQ( points(Y,i), y[i] );
    }
    
}

TEST_F(GeometryTest, PlungingBody) {
    // Define a flat plate with 10 points
    int nPoints = 10;
    RigidBody body;
    body.addLine_n( 0, 0, 1, 0, nPoints );
    EXPECT_EQ( true, body.isStationary() );

    // Add it to a geometry
    Geometry* geom = new Geometry;
    geom->addBody(body);
    geom->moveBodies(0);
    EXPECT_EQ( true, geom->isStationary() );
    
    // Check velocities are zero
    BoundaryVector vel = geom->getVelocities();
    EXPECT_EQ( nPoints, vel.getNumPoints() );
    for (int i=0; i < nPoints; ++i) {
        EXPECT_DOUBLE_EQ( 0, vel(X,i) );
        EXPECT_DOUBLE_EQ( 0, vel(Y,i) );
    }
    delete geom;

    // Define a plunging motion y(t) = A sin(2pi w t)
    double amp = 2;
    double freq = 3;
    PitchPlunge motion( 0, 0, amp, freq );

    // Add it to the body
    body.setMotion( motion );
    EXPECT_EQ( false, body.isStationary() );

    // Add the body to the geometry
    geom = new Geometry;
    geom->addBody(body);
    EXPECT_EQ( false, geom->isStationary() );

    // Check velocities are ( 0, amp * 2pi * freq) at time 0
    geom->moveBodies(0);
    vel = geom->getVelocities();
    EXPECT_EQ( nPoints, vel.getNumPoints() );
    for (int i=0; i < nPoints; ++i) {
        EXPECT_DOUBLE_EQ( 0, vel(X,i) );
        EXPECT_DOUBLE_EQ( amp * 2 * PI * freq, vel(Y,i) );
    }
    
    // Check velocities are ( 0, 0 ) at time T/4
    double period = 1. / freq;
    geom->moveBodies( period / 4. );
    vel = geom->getVelocities();
    EXPECT_EQ( nPoints, vel.getNumPoints() );
    for (int i=0; i < nPoints; ++i) {
        EXPECT_NEAR( 0, vel(X,i), 1e-14 );
        EXPECT_NEAR( 0, vel(Y,i), 1e-14 );
    }

    delete geom;
}

}  // namespace

