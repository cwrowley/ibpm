#include "OutputProbes.h"
#include <gtest/gtest.h>

using namespace ibpm;

namespace {

double tol = 1e-14;
    
class OutputProbesTest : public testing::Test {
protected:
    OutputProbesTest( ) : 
        _nx(8),
        _ny(12),
        _ngrid(3),
        _xLength(2),
        _yLength(_xLength * _ny/_nx),
        _xOffset(-1),
        _yOffset(-3),
        _dx(_xLength/_nx),
        _filename("testOutputProbes.probe"), 
        _grid(_nx, _ny, _ngrid, _xLength, _xOffset, _yOffset),
		_flagInitialization(false),
		_probe(_filename, _grid) {			
		}
	
    virtual ~OutputProbesTest() { }

    int _nx;
    int _ny;
    int _ngrid;
    double _xLength;
    double _yLength;
    double _xOffset;
    double _yOffset;
    double _dx;
	string _filename;
	Grid _grid;
	int _lev;
	int _numProbes;
	int _dimen;
	bool _flagInitialization;
	Array::Array2<int> _probePositions;
	OutputProbes _probe;
};

TEST_F(OutputProbesTest, TestgetProbeInfoFunctions){
// test getNumProbes, getProbeIndexX, getProbeIndexY, addProbebyIndex, addProbebyPosition
	int i = 4; 	int j = 3;  // probe 1
	int ii = 6; int jj = 8; // probe 2
	int lev = 0;
	OutputProbes probe(_filename, _grid);
	EXPECT_DOUBLE_EQ( 0, probe.getNumProbes() );		
	probe.addProbeByIndex( i , j ); 
	EXPECT_DOUBLE_EQ( 1, probe.getNumProbes() );
	probe.addProbeByIndex( ii, jj ); 
	EXPECT_DOUBLE_EQ( 2, probe.getNumProbes() );
	EXPECT_DOUBLE_EQ( i, probe.getProbeIndexX(1) );
	EXPECT_DOUBLE_EQ( j, probe.getProbeIndexY(1) );
	EXPECT_DOUBLE_EQ( ii, probe.getProbeIndexX(2) );
	EXPECT_DOUBLE_EQ( jj, probe.getProbeIndexY(2) );
	EXPECT_DOUBLE_EQ( _grid.getXEdge(lev, ii), probe.getProbeCoordX(2) );
	EXPECT_DOUBLE_EQ( _grid.getYEdge(lev, jj), probe.getProbeCoordY(2) );
	probe.addProbeByPosition( _grid.getXEdge(lev, ii), _grid.getYEdge(lev, jj) ); 
	EXPECT_DOUBLE_EQ( 3, probe.getNumProbes() );
	EXPECT_DOUBLE_EQ( ii, probe.getProbeIndexX(3) );
	EXPECT_DOUBLE_EQ( jj, probe.getProbeIndexY(3) );
}

}  // namespace
