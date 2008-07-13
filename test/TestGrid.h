#include <cxxtest/TestSuite.h>
#include "Grid.h"

class TestGrid : public CxxTest::TestSuite {
public:
	void setUp() {
		int nx = 10;
		int ny = 20;
		double length = 2;
		double xOffset = -1;
		double yOffset = -3;
		_grid = new Grid(nx, ny, length, xOffset, yOffset);
	}

	void tearDown() {
		delete _grid;
	}

	void testDx() {
		TS_ASSERT_DELTA(_grid->getDx(), 0.2, _delta);
	}

	void testXCenter() {
		TS_ASSERT_DELTA(_grid->getXCenter(0),-0.9, _delta);
		TS_ASSERT_DELTA(_grid->getXCenter(9), 0.9, _delta);	
		TS_ASSERT_DELTA(_grid->getXCenter(5), 0.1, _delta);
	}

	void testYCenter() {
		TS_ASSERT_DELTA(_grid->getYCenter(0), -2.9, _delta);
		TS_ASSERT_DELTA(_grid->getYCenter(10),-0.9, _delta);	
		TS_ASSERT_DELTA(_grid->getYCenter(19), 0.9, _delta);
	}

	void testXEdge() {
		TS_ASSERT_DELTA(_grid->getXEdge(0),-1, _delta);
		TS_ASSERT_DELTA(_grid->getXEdge(5), 0, _delta);
		TS_ASSERT_DELTA(_grid->getXEdge(10),1, _delta);		
	}

	void testYEdge() {
		TS_ASSERT_DELTA(_grid->getYEdge(0), -3, _delta);
		TS_ASSERT_DELTA(_grid->getYEdge(10),-1, _delta);
		TS_ASSERT_DELTA(_grid->getYEdge(20), 1, _delta);		
	}
    // void testBoundaries() {
    //  TS_ASSERT_THROWS_ANYTHING(_grid->getXCenter(-1));
    //  TS_ASSERT_THROWS_ANYTHING(_grid->getYCenter(-1));
    //  TS_ASSERT_THROWS_ANYTHING(_grid->getXCenter(10));
    //  TS_ASSERT_THROWS_ANYTHING(_grid->getYCenter(20));
    //  TS_ASSERT_THROWS_ANYTHING(_grid->getXEdge(-1));
    //  TS_ASSERT_THROWS_ANYTHING(_grid->getYEdge(-1));
    //  TS_ASSERT_THROWS_ANYTHING(_grid->getXEdge(11));
    //  TS_ASSERT_THROWS_ANYTHING(_grid->getYEdge(21));
    // }

private:
	Grid* _grid;
	static double _delta;  // error tolerance for floating-point tests
};

double TestGrid::_delta = 1e-10;