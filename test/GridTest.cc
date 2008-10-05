#include "Grid.h"
#include <gtest/gtest.h>

using namespace ibpm;

namespace {

double tol = 1e-14;
    
class GridTest : public testing::Test {
protected:
    GridTest() : 
        _nx(8),
        _ny(12),
        _ngrid(3),
        _xLength(2),
        _yLength(_xLength * _ny/_nx),
        _xOffset(-1),
        _yOffset(-3),
        _dx(_xLength/_nx),
        _grid(_nx, _ny, _ngrid, _xLength, _xOffset, _yOffset) {}

    virtual ~GridTest() {}

    int _nx;
    int _ny;
    int _ngrid;
    double _xLength;
    double _yLength;
    double _xOffset;
    double _yOffset;
    double _dx;
    Grid _grid;
};

TEST_F(GridTest, TestDx) {
    EXPECT_DOUBLE_EQ( _dx, _grid.Dx() );
}

TEST_F(GridTest, TestXCenter) {
     EXPECT_DOUBLE_EQ( _xOffset + _dx/2, _grid.getXCenter(0,0) );
     EXPECT_DOUBLE_EQ( _xOffset + _xLength - _dx/2, _grid.getXCenter(0,_nx-1)); 
     EXPECT_NEAR( _xOffset + _xLength/2 + _dx/2, _grid.getXCenter(0,_nx/2), tol);
}

TEST_F(GridTest, TestYCenter) {
    EXPECT_DOUBLE_EQ( _yOffset + _dx/2, _grid.getYCenter(0,0) );
    EXPECT_DOUBLE_EQ( _yOffset + _yLength - _dx/2, _grid.getYCenter(0,_ny-1)); 
    EXPECT_NEAR( _yOffset + _yLength/2 + _dx/2, _grid.getYCenter(0,_ny/2), tol);
}
    
TEST_F(GridTest, TestXBoundaries) {
    EXPECT_DOUBLE_EQ( _xOffset, _grid.getXEdge(0,0) );
    EXPECT_DOUBLE_EQ( _xOffset - _xLength/2, _grid.getXEdge(1,0) );
    EXPECT_DOUBLE_EQ( _xOffset - 3*_xLength/2, _grid.getXEdge(2,0));
    EXPECT_DOUBLE_EQ( _xLength+_xOffset, _grid.getXEdge(0,_nx) );
    EXPECT_DOUBLE_EQ( 3*_xLength/2 + _xOffset, _grid.getXEdge(1,_nx) );
    EXPECT_DOUBLE_EQ( 5*_xLength/2 + _xOffset, _grid.getXEdge(2,_nx) );
}                         
    
TEST_F(GridTest, TestYBoundaries) {
    EXPECT_DOUBLE_EQ( _yOffset, _grid.getYEdge(0,0) );
    EXPECT_DOUBLE_EQ( _yOffset - _yLength/2, _grid.getYEdge(1,0) );
    EXPECT_DOUBLE_EQ( _yOffset - 3*_yLength/2, _grid.getYEdge(2,0));
    EXPECT_DOUBLE_EQ( _yLength+_yOffset, _grid.getYEdge(0,_ny) );
    EXPECT_DOUBLE_EQ( 3*_yLength/2 + _yOffset, _grid.getYEdge(1,_ny) );
    EXPECT_DOUBLE_EQ( 5*_yLength/2 + _yOffset, _grid.getYEdge(2,_ny) );
}                         


//TEST_F(GridTest, testXEdge) {
//    EXPECT_DOUBLE_EQ(_grid.getXEdge(0,0),-1);
//    EXPECT_DOUBLE_EQ(_grid.getXEdge(0,5), 0);
//    EXPECT_DOUBLE_EQ(_grid.getXEdge(0,10),1);
//    EXPECT_DOUBLE_EQ(_grid.getXEdge(1,0), _xOffset - _xLength/4);
//    EXPECT_DOUBLE_EQ(_grid.getXEdge(2,0), _xOffset - _xLength/4 - _xLength/2);
//}
//
//TEST_F(GridTest, testYEdge) {
//    EXPECT_DOUBLE_EQ(_grid.getYEdge(0,0), -3);
//    EXPECT_DOUBLE_EQ(_grid.getYEdge(0,10),-1);
//    EXPECT_DOUBLE_EQ(_grid.getYEdge(0,20), 1);        
//}

}  // namespace

