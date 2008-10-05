#include "State.h"
#include <unistd.h>
#include <gtest/gtest.h>

using namespace ibpm;

namespace {

#define EXPECT_SCALAR_EQ(a,b)                   \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for (int i=1; i<_nx; ++i) {             \
            for (int j=1; j<_ny; ++j) {         \
                EXPECT_DOUBLE_EQ((a), (b));     \
            }                                   \
        }                                       \
    }
    
#define EXPECT_FLUX_X_EQ(a,b)                   \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for ( int i=0; i<_nx+1; ++i ) {         \
            for ( int j=0; j<_ny; ++j ) {       \
                EXPECT_DOUBLE_EQ( (a), (b) );   \
            }                                   \
        }                                       \
    }

#define EXPECT_FLUX_Y_EQ(a,b)                   \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for ( int i=0; i<_nx; ++i ) {           \
            for ( int j=0; j<_ny+1; ++j ) {     \
                EXPECT_DOUBLE_EQ( (a), (b) );   \
            }                                   \
        }                                       \
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
        _nx(4),
        _ny(8),
        _ngrid(3),
        _numPoints(2),
        _grid( _nx, _ny, _ngrid, 2, -1, -2 ),
        _x( _grid, _numPoints ) {
        
        // Initialize a State
        _x.q = 1.;
        _x.omega = 2.;
        _x.f = 3.;
        _x.timestep = 4;
        _x.time = 0.5;
    }

    // data
    int _nx;
    int _ny;
    int _ngrid;
    int _numPoints;
    Grid _grid;
    State _x;
};

TEST_F( StateTest, Save ) {
    bool success = _x.save("state_test");
    EXPECT_EQ( true, success );
    
    State y( _grid, _numPoints );
    y.load( "state_test" );
    EXPECT_SCALAR_EQ( _x.omega(lev,i,j), y.omega(lev,i,j) );
    EXPECT_FLUX_X_EQ( _x.q(lev,X,i,j),   y.q(lev,X,i,j)   );
    EXPECT_FLUX_Y_EQ( _x.q(lev,Y,i,j),   y.q(lev,Y,i,j)   );
    EXPECT_BV_EQ(     _x.f(dir,i),   y.f(dir,i)   );
    EXPECT_DOUBLE_EQ( _x.time,       y.time       );
    EXPECT_EQ(        _x.timestep,   y.timestep   );
    
    unlink("state_test");
}

TEST_F( StateTest, LoadBadGeometry ) {
    bool status = _x.save("state_test");
    ASSERT_EQ( true, status );

    // Create another state with a different Geometry (no points)
    State y( _grid, 0 );
    status = y.load( "state_test" );
    EXPECT_EQ( false, status );

    unlink("state_test");
}

TEST_F( StateTest, LoadBadGrid ) {
    bool status = _x.save("state_test");
    ASSERT_EQ( true, status );

    // Create another state with a different Grid
    Grid newGrid( 4, 5, 1, 2, -1, -2 );
    State y( newGrid, _numPoints );
    status = y.load( "state_test" );
    EXPECT_EQ( false, status );
    
    unlink("state_test");
}

} // namespace
