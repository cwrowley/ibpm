#include "Grid.h"
#include <gtest/gtest.h>
#include <assert.h>

using namespace ibpm;

namespace {

const double tol = 1e-14;

// Use value parameterized tests to test the grid shifting methods
// The shift values must be global for the parameterized tests 
// _shiftVal contains the *number of grid points* by which grid level 1 is shifted 
// _shiftVal cannot take magnitudes greater than _nxShort/4 (defined below)
const int _shiftVal[] = {-2., -1., 0., 1., 2.};	
	
// Each test will be performed for a square, fat, and tall grid
// Note that each grid is set up so that _dx is the same for each case
class GridTest : public testing::TestWithParam<int> {
protected:
	GridTest() : 
    _nxLong(12),
    _nxShort(8),
	_ngrid(3),
	_xLengthLong(2),
	_xLengthShort(_xLengthLong * _nxShort/_nxLong),
	_xOffset(-1),
	_yOffset(-3),
    _dx(_xLengthLong / _nxLong),
    _allGrids(3) {
        // Square grid
        _allGrids[0] = Grid(_nxLong, _nxLong, _ngrid, _xLengthLong, _xOffset, _yOffset);
        // Fat grid
        _allGrids[1] = Grid(_nxLong, _nxShort, _ngrid, _xLengthLong, _xOffset, _yOffset);
        // Tall grid
        _allGrids[2] = Grid(_nxShort, _nxLong, _ngrid, _xLengthShort, _xOffset, _yOffset);
    }
	
	virtual ~GridTest() {}
	
	int _nxLong;
    int _nxShort;
    int _ngrid;
	double _xLengthLong;
	double _xLengthShort;
	double _xOffset;
	double _yOffset;
    double _dx;
    vector< Grid > _allGrids;
    vector< Grid >::iterator _grid;  
};
		
TEST_F(GridTest, TestDx) {
    for ( _grid = _allGrids.begin(); _grid != _allGrids.end(); _grid++) {
        EXPECT_DOUBLE_EQ( _dx, _grid->Dx() );
        EXPECT_DOUBLE_EQ( _dx, _grid->Dx(0) );
        EXPECT_DOUBLE_EQ( _dx * 2, _grid->Dx(1) );
        EXPECT_DOUBLE_EQ( _dx * 4, _grid->Dx(2) );  
   }
}

TEST_F(GridTest, TestXCenter) {
	for ( _grid = _allGrids.begin(); _grid != _allGrids.end(); _grid++) {
		int _nx = _grid->Nx();
		double _xLength = _nx * _dx;
		EXPECT_DOUBLE_EQ( _xOffset + _dx/2, _grid->getXCenter(0,0) );
		EXPECT_DOUBLE_EQ( _xOffset + _xLength - _dx/2, _grid->getXCenter(0,_nx-1)); 
		EXPECT_NEAR( _xOffset + _xLength/2 + _dx/2, _grid->getXCenter(0,_nx/2), tol);
	}
}

TEST_F(GridTest, TestYCenter) {
	for ( _grid = _allGrids.begin(); _grid != _allGrids.end(); _grid++) {
		int _ny = _grid->Ny();
		double _yLength = _ny * _dx;
		EXPECT_DOUBLE_EQ( _yOffset + _dx/2, _grid->getYCenter(0,0) );
		EXPECT_DOUBLE_EQ( _yOffset + _yLength - _dx/2, _grid->getYCenter(0,_ny-1)); 
		EXPECT_NEAR( _yOffset + _yLength/2 + _dx/2, _grid->getYCenter(0,_ny/2), tol);
	}
}                      
 
// This tests the grid boundaries at two grid levels for x-shifts
TEST_P(GridTest, TestXShift) {
	for ( _grid = _allGrids.begin(); _grid != _allGrids.end(); _grid++) {
		int _nx = _grid->Nx();
		double _xLength = (double) _nx * _dx;
		
		int nShift = GetParam();
		double xShiftVal = (double) nShift*4/_nx;
		
		// Check x-shift
		_grid->setXShift( xShiftVal );

		//Check left boundaries
		EXPECT_NEAR( _xOffset, _grid->getXEdge(0,0), tol );
		EXPECT_NEAR( _xOffset - (_nx/4-nShift)*_grid->Dx(1), _grid->getXEdge(1,0), tol );
		EXPECT_NEAR( _grid->getXEdge(1,0) - (_nx/4-nShift)*_grid->Dx(2), _grid->getXEdge(2,0), tol );
		
		// Check right boundaries
		EXPECT_NEAR( _xLength+_xOffset, _grid->getXEdge(0,_nx), tol );
		EXPECT_NEAR( _xOffset + (3*_nx/4+nShift)*_grid->Dx(1), _grid->getXEdge(1,_nx), tol );
		EXPECT_NEAR( _grid->getXEdge(1,0) + (3*_nx/4+nShift)*_grid->Dx(2), _grid->getXEdge(2,_nx), tol );
	}
}
	
// Test grid boundaries for y-shifts
TEST_P(GridTest, TestYShift) {
	for ( _grid = _allGrids.begin(); _grid != _allGrids.end(); _grid++) {
		int _ny = _grid->Ny();
		double _yLength = (double) _ny * _dx;
		
		int nShift = GetParam();
		double yShiftVal = (double) nShift*4/_ny;
		
		// Check y-shift
		_grid->setYShift( yShiftVal );
		
		//Check bottom boundaries
		EXPECT_NEAR( _yOffset, _grid->getYEdge(0,0), tol );
		EXPECT_NEAR( _yOffset - (_ny/4-nShift)*_grid->Dx(1), _grid->getYEdge(1,0), tol );
		EXPECT_NEAR( _grid->getYEdge(1,0) - (_ny/4-nShift)*_grid->Dx(2), _grid->getYEdge(2,0), tol );
		
		// Check top boundaries
		EXPECT_NEAR( _yLength+_yOffset, _grid->getYEdge(0,_ny), tol );
		EXPECT_NEAR( _yOffset + (3*_ny/4+nShift)*_grid->Dx(1), _grid->getYEdge(1,_ny), tol );
		EXPECT_NEAR( _grid->getYEdge(1,0) + (3*_ny/4+nShift)*_grid->Dx(2), _grid->getYEdge(2,_ny), tol );
	}
}
	
INSTANTIATE_TEST_CASE_P(
	shiftTests, GridTest, ::testing::ValuesIn(_shiftVal) 
);

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

