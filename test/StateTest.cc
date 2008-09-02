#include "State.h"
#include <unistd.h>
#include <gtest/gtest.h>

using namespace ibpm;

namespace {

#define EXPECT_SCALAR_EQ(a,b)               \
    for (int i=0; i<_nx+1; ++i) {           \
        for (int j=0; j<_ny+1; ++j) {       \
            EXPECT_DOUBLE_EQ((a), (b));     \
        }                                   \
    }

#define EXPECT_FLUX_X_EQ(a,b)             \
    for ( int i=0; i<_nx+1; ++i ) {       \
        for ( int j=0; j<_ny; ++j ) {     \
            EXPECT_DOUBLE_EQ( (a), (b) ); \
        }                                 \
    }

#define EXPECT_FLUX_Y_EQ(a,b)             \
    for ( int i=0; i<_nx; ++i ) {         \
        for ( int j=0; j<_ny+1; ++j ) {   \
            EXPECT_DOUBLE_EQ( (a), (b) ); \
        }                                 \
    }

#define EXPECT_BV_EQ(a,b)                       \
    for ( Direction dir=X; dir < XY; ++dir ) {  \
        for ( int i=0; i<_numPoints; ++i ) {    \
            EXPECT_DOUBLE_EQ( (a), (b) );       \
        }                                       \
    }

class StateTest : public testing::Test {
protected:
    StateTest() :
        _nx(3),
        _ny(5),
        _numPoints(1),
        _grid( _nx, _ny, 2, -1, -2 ) {
        
        // Create a geometry with one point
        RigidBody body;
        body.addPoint(0,0);
        _geom.addBody(body);
        _x = new State( _grid, _geom );

        // Initialize a State
        _x->q = 1.;
        _x->gamma = 2.;
        _x->f = 3.;
        _x->timestep = 4;
        _x->time = 0.5;
    }

    // data
    int _nx;
    int _ny;
    int _numPoints;
    Grid _grid;
    Geometry _geom;
    State* _x;
};

TEST_F( StateTest, Save ) {
    bool success = _x->save("state_test");
    EXPECT_EQ( true, success );
    
    State y(_grid, _geom);
    y.load( "state_test" );
    EXPECT_SCALAR_EQ( _x->gamma(i,j), y.gamma(i,j) );
    EXPECT_FLUX_X_EQ( _x->q(X,i,j),   y.q(X,i,j)   );
    EXPECT_FLUX_Y_EQ( _x->q(Y,i,j),   y.q(Y,i,j)   );
    EXPECT_BV_EQ(     _x->f(dir,i),   y.f(dir,i)   );
    EXPECT_DOUBLE_EQ( _x->time,       y.time       );
    EXPECT_EQ(        _x->timestep,   y.timestep   );
    EXPECT_DOUBLE_EQ( _x->q(X,0,0),   y.q(X,0,0)   );
    
    unlink("state_test");
}

TEST_F( StateTest, LoadBadGeometry ) {
    bool status = _x->save("state_test");
    ASSERT_EQ( true, status );

    // Create another state with a different Geometry
    Geometry emptyGeom;
    State y( _grid, emptyGeom );
    status = y.load( "state_test" );
    EXPECT_EQ( false, status );

    unlink("state_test");
}

TEST_F( StateTest, LoadBadGrid ) {
    bool status = _x->save("state_test");
    ASSERT_EQ( true, status );

    // Create another state with a different Grid
    Grid newGrid( 4, 5, 2, -1, -2 );
    State y( newGrid, _geom );
    status = y.load( "state_test" );
    EXPECT_EQ( false, status );
    
    unlink("state_test");
}

} // namespace