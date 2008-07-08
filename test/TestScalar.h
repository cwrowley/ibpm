#include <cxxtest/TestSuite.h>
#include "Grid.h"
#include "Scalar.h"

class TestScalar : public CxxTest::TestSuite {
public:
	void setUp() {
		// define a Grid
		_nx = 3;
		_ny = 5;
		double length = 2;
		double xOffset = -1;
		double yOffset = -3;
		_grid = new Grid( _nx, _ny, length, xOffset, yOffset );

		_f = new Scalar( *_grid );
		_g = new Scalar( *_grid );
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny; ++j) {
				(*_f)(i,j) = f(i,j);
				(*_g)(i,j) = g(i,j);
			}
		}
	}

	void testAssignment() {
		int i = _nx/2;
		int j = _ny/2;
		(*_f)(i,j) = 73;
		TS_ASSERT_DELTA((*_f)(i,j), 73, _delta);
	}

	void testCopyConstructor() {
		Scalar h = *_f;
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny; ++j) {
				TS_ASSERT_DELTA(h(i,j), f(i,j), _delta);
			}
		}
	}

	void testCopyFromDouble() {
		double a = 7;
		*_f = a;
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny; ++j) {
				TS_ASSERT_DELTA( (*_f)(i,j), 7, _delta);
			}
		}		
	}

	void testElements() {
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny; ++j) {
				TS_ASSERT_DELTA((*_f)(i,j), f(i,j), _delta);
			}
		}
	}

	void testPlusEqualsScalar() {
		*_g += *_f;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( (*_g)(i,j), f(i,j) + g(i,j), _delta );
			}
		}
	}

	void testPlusEqualsDouble() {
		*_g += 6;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( (*_g)(i,j), g(i,j)+6, _delta );
			}
		}
	}

	void testScalarPlusScalar() {
		Scalar h = (*_f) + (*_g);
		
		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), f(i,j) + g(i,j), _delta );
			}
		}
	}

	void testScalarPlusDouble() {
		Scalar h = (*_f) + 4;
		
		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), f(i,j) + 4, _delta );
			}
		}
	}

	void testDoublePlusScalar() {
		Scalar h = 4 + (*_f);
		
		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), f(i,j) + 4, _delta );
			}
		}
	}

	void testMinusEquals() {
		*_g -= *_f;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( (*_g)(i,j), g(i,j) - f(i,j), _delta );
			}
		}
	}

	void testScalarMinusScalar() {
		Scalar h = *_g - *_f;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), g(i,j) - f(i,j), _delta );
			}
		}
	}
	
	void testScalarMinusDouble() {
		Scalar h = *_g - 9;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), g(i,j) - 9, _delta );
			}
		}
	}
	
	void testDoubleMinusScalar() {
		Scalar h = 9 - *_g;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), 9 - g(i,j), _delta );
			}
		}
	}
	
	void testTimesEquals() {
		*_g *= *_f;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( (*_g)(i,j), g(i,j) * f(i,j), _delta );
			}
		}
	}

	void testTimesEqualsDouble() {
		*_g *= 3;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( (*_g)(i,j), g(i,j) * 3, _delta );
			}
		}
	}

	void testScalarTimesScalar() {
		Scalar h = (*_g) * (*_f);

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), g(i,j) * f(i,j), _delta );
			}
		}
	}
	
	void testScalarTimesDouble() {
		Scalar h = (*_g) * 11;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), g(i,j) * 11, _delta );
			}
		}
	}
	
	void testDoubleTimesScalar() {
		Scalar h = 11 * (*_g);

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), g(i,j) * 11, _delta );
			}
		}
	}
	
	void testDivEquals() {
		*_g /= *_f;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( (*_g)(i,j), g(i,j) / f(i,j), _delta );
			}
		}
	}

	void testDivEqualsDouble() {
		*_g /= 3;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( (*_g)(i,j), g(i,j) / 3, _delta );
			}
		}
	}

	void testScalarDivScalar() {
		Scalar h = (*_g) / (*_f);

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), g(i,j) / f(i,j), _delta );
			}
		}
	}
	
	void testScalarDivDouble() {
		Scalar h = (*_g) / 11;

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), g(i,j) / 11, _delta );
			}
		}
	}
	
	void testDoubleDivScalar() {
		Scalar h = 11 / (*_g);

		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), 11 / g(i,j), _delta );
			}
		}
	}
	
	void testUnaryMinus() {
		Scalar h = -(*_f);
		
		for ( int i=0; i<_nx; ++i ) {
			for ( int j=0; j<_ny; ++j ) {
				TS_ASSERT_DELTA( h(i,j), -f(i,j), _delta );
			}
		}
	}
	
	void tearDown() {
		delete _f;
	}

	// functions to test
	double f(int i, int j) {
		return i*_nx + j + 1;
	}

	double g(int i, int j) {
		return 3*i + j;
	}

private:
	int _nx;
	int _ny;
	Scalar* _f;
	Scalar* _g;
	Grid* _grid;
	static double _delta;
};

double TestScalar::_delta = 1e-10;
