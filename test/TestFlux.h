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