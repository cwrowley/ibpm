#include "BC.h"
#include <gtest/gtest.h>

using namespace ibpm;

class BCTest : public testing::Test {
protected:
    BCTest() :
        _nx(8),
        _ny(4),
        _bc( _nx, _ny )
    {}
    
    ~BCTest() {}

    // data
    int _nx;
    int _ny;
    BC _bc;
};

TEST_F( BCTest, InitializeToZero ) {
    for ( int i=0; i<=_nx; ++i ) {
        EXPECT_DOUBLE_EQ( 0., _bc.bottom(i) );
        EXPECT_DOUBLE_EQ( 0., _bc.top(i) );
    }
    for ( int j=0; j <= _ny; ++j ) {
        EXPECT_DOUBLE_EQ( 0., _bc.left(j) );
        EXPECT_DOUBLE_EQ( 0., _bc.right(j) );
    }
}

TEST_F( BCTest, CheckAssignment ) {
    // left boundary
    _bc.left(0) = 3.;
    EXPECT_DOUBLE_EQ( 3., _bc.left(0) );
    _bc.left(_ny/2) = 4.;
    EXPECT_DOUBLE_EQ( 4., _bc.left(_ny/2) );
    _bc.left(_ny) = 5.;
    EXPECT_DOUBLE_EQ( 5., _bc.left(_ny) );
    
    // right boundary
    _bc.right(0) = 6.;
    EXPECT_DOUBLE_EQ( 6., _bc.right(0) );
    _bc.right(_ny/2) = 7.;
    EXPECT_DOUBLE_EQ( 7., _bc.right(_ny/2) );
    _bc.right(_ny) = 8.;
    EXPECT_DOUBLE_EQ( 8., _bc.right(_ny) );

    // top boundary
    _bc.top(0) = 9.;
    EXPECT_DOUBLE_EQ( 9., _bc.top(0) );
    _bc.top(_nx/2) = 10.;
    EXPECT_DOUBLE_EQ( 10., _bc.top(_nx/2) );
    _bc.top(_nx) = 11.;
    EXPECT_DOUBLE_EQ( 11., _bc.top(_nx) );
    
    // bottom boundary
    _bc.bottom(0) = 12.;
    EXPECT_DOUBLE_EQ( 12., _bc.bottom(0) );
    _bc.bottom(_nx/2) = 13.;
    EXPECT_DOUBLE_EQ( 13., _bc.bottom(_nx/2) );
    _bc.bottom(_nx) = 14.;
    EXPECT_DOUBLE_EQ( 14., _bc.bottom(_nx) );
    
}

TEST_F( BCTest, CheckCorners ) {
    // left(0) should equal bottom(0)
    _bc.left(0) = 3.;
    EXPECT_DOUBLE_EQ( 3., _bc.bottom(0) );
    _bc.bottom(0) = 4.;
    EXPECT_DOUBLE_EQ( 4., _bc.left(0) );

    // top(0) should equal left(ny)
    _bc.top(0) = 6.;
    EXPECT_DOUBLE_EQ( 6., _bc.left(_ny) );
    _bc.left(_ny) = 7.;
    EXPECT_DOUBLE_EQ( 7., _bc.top(0) );
    
    // top(nx) should equal right(ny)
    _bc.right(_ny) = 8;
    EXPECT_DOUBLE_EQ( 8., _bc.top(_nx) );
    _bc.top(_nx) = 9.;
    EXPECT_DOUBLE_EQ( 9., _bc.right(_ny) );

    // bottom(nx) should equal right(0)
    _bc.bottom(_nx) = 10.;
    EXPECT_DOUBLE_EQ( 10., _bc.right(0) );
    _bc.right(0) = 11.;
    EXPECT_DOUBLE_EQ( 11., _bc.bottom(_nx) );
}

TEST_F( BCTest, BCEqualsDouble ) {
    const double val = 7.;
    _bc = val;
    for (int i=0; i<=_nx; ++i) {
        EXPECT_DOUBLE_EQ( val, _bc.bottom(i) );
        EXPECT_DOUBLE_EQ( val, _bc.top(i) );
    }
    for (int j=0; j<=_ny; ++j) {
        EXPECT_DOUBLE_EQ( val, _bc.left(j) );
        EXPECT_DOUBLE_EQ( val, _bc.right(j) );
    }    
}

