#include <cxxtest/TestSuite.h>
#include "Grid.h"
#include "Flux.h"

class TestFlux : public CxxTest::TestSuite {
public:
	void setUp() {
		// define a Grid
		int nx = 10;
		int ny = 20;
		double length = 2;
		double xOffset = -1;
		double yOffset = -3;
		Grid grid( nx, ny, length, xOffset, yOffset );

		_f = new Flux( grid );
		
		// Is grid dimension equal to flux dimension?
		TS_ASSERT_DELTA(grid.getNx(),_f->getNx(),_delta);
		TS_ASSERT_DELTA(grid.getNy(),_f->getNy(),_delta);
   		// Is grid dimension correct?
		TS_ASSERT_DELTA(grid.getNx(),nx,_delta);
		TS_ASSERT_DELTA(grid.getNy(),ny,_delta);
	}

	void tearDown() {
		delete _f;
	}

	void testTrue() {
		TS_ASSERT_DELTA(1, 1, _delta);
	}
	
	// WORK HERE: add some real tests

private:
	Flux* _f;
	static double _delta;
};

double TestFlux::_delta = 1e-10;
