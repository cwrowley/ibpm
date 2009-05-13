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
    EXPECT_DOUBLE_EQ( _dx, _grid.Dx(0) );
    EXPECT_DOUBLE_EQ( _dx * 2, _grid.Dx(1) );
    EXPECT_DOUBLE_EQ( _dx * 4, _grid.Dx(2) );    
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

TEST_F(GridTest, TestXIndex) {
	double i = 5;
	double drift1 = (_dx / 2 ) * 0.94;
	double x1 = i * _dx + drift1;
	double x2 = (i + 1) * _dx - drift1;
		
	EXPECT_DOUBLE_EQ( 0, _grid.getXGridIndex(_xOffset) );
	EXPECT_DOUBLE_EQ( _nx, _grid.getXGridIndex(_xOffset + _xLength) );
	EXPECT_DOUBLE_EQ( 1, _grid.getXGridIndex(_xOffset + _dx) );	
	EXPECT_DOUBLE_EQ( _nx - 1, _grid.getXGridIndex(_xOffset + _xLength - _dx) );
	EXPECT_DOUBLE_EQ( i , _grid.getXGridIndex(_xOffset + x1) );
	EXPECT_DOUBLE_EQ( i + 1, _grid.getXGridIndex(_xOffset + x2) );
}

TEST_F(GridTest, TestYIndex) {
	double j = 5;
	double drift1 = (_dx / 2 ) * 0.94;
	double y1 = j * _dx + drift1;
	double y2 = (j + 1) * _dx - drift1;
		
	EXPECT_DOUBLE_EQ( 0, _grid.getYGridIndex(_yOffset) );
	EXPECT_DOUBLE_EQ( _ny, _grid.getYGridIndex(_yOffset + _yLength) );
	EXPECT_DOUBLE_EQ( 1, _grid.getYGridIndex(_yOffset + _dx) );	
	EXPECT_DOUBLE_EQ( _ny - 1, _grid.getYGridIndex(_yOffset + _yLength - _dx) );
	EXPECT_DOUBLE_EQ( j , _grid.getYGridIndex(_yOffset + y1) );
	EXPECT_DOUBLE_EQ( j + 1, _grid.getYGridIndex(_yOffset + y2) );
}

}  // namespace

