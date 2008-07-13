#include <cxxtest/TestSuite.h>
#include "Grid.h"
#include "Scalar.h"

#define ASSERT_ALL_EQUAL(a,b)                        \
	for (int i=0; i<_nx; ++i) {                      \
		for (int j=0; j<_ny; ++j) {                  \
			TS_ASSERT_DELTA((a), (b), _delta);       \
		}                                            \
	}

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
		ASSERT_ALL_EQUAL( h(i,j), f(i,j) );
	}

	void testCopyFromDouble() {
		double a = 7;
		*_f = a;
		ASSERT_ALL_EQUAL( (*_f)(i,j), a );
	}

	void testElements() {
		ASSERT_ALL_EQUAL( (*_f)(i,j), f(i,j) );
	}

	void testPlusEqualsScalar() {
		*_g += *_f;
		ASSERT_ALL_EQUAL( (*_g)(i,j), f(i,j) + g(i,j) );
	}

	void testPlusEqualsDouble() {
		*_g += 6;
		ASSERT_ALL_EQUAL( (*_g)(i,j), g(i,j) + 6 );
	}

	void testScalarPlusScalar() {
		Scalar h = (*_f) + (*_g);
		ASSERT_ALL_EQUAL( h(i,j), f(i,j) + g(i,j) );
	}

	void testScalarPlusDouble() {
		Scalar h = (*_f) + 4;
		ASSERT_ALL_EQUAL( h(i,j), f(i,j) + 4 );
	}

	void testDoublePlusScalar() {
		Scalar h = 4 + (*_f);
		ASSERT_ALL_EQUAL( h(i,j), f(i,j) + 4 );
	}

	void testMinusEquals() {
		*_g -= *_f;
		ASSERT_ALL_EQUAL( (*_g)(i,j), g(i,j) - f(i,j) );
	}

	void testScalarMinusScalar() {
		Scalar h = *_g - *_f;
		ASSERT_ALL_EQUAL( h(i,j), g(i,j) - f(i,j) );
	}
	
	void testScalarMinusDouble() {
		Scalar h = *_g - 9;
		ASSERT_ALL_EQUAL( h(i,j), g(i,j) - 9 );
	}
	
	void testDoubleMinusScalar() {
		Scalar h = 9 - *_g;
		ASSERT_ALL_EQUAL( h(i,j), 9 - g(i,j) );
	}
	
	void testTimesEquals() {
		*_g *= *_f;
		ASSERT_ALL_EQUAL( (*_g)(i,j), g(i,j) * f(i,j) );
	}

	void testTimesEqualsDouble() {
		*_g *= 3;
		ASSERT_ALL_EQUAL( (*_g)(i,j), g(i,j) * 3 );
	}

	void testScalarTimesScalar() {
		Scalar h = (*_g) * (*_f);
		ASSERT_ALL_EQUAL( h(i,j), g(i,j) * f(i,j) );
	}
	
	void testScalarTimesDouble() {
		Scalar h = (*_g) * 11;
		ASSERT_ALL_EQUAL( h(i,j), g(i,j) * 11 );
	}
	
	void testDoubleTimesScalar() {
		Scalar h = 11 * (*_g);
		ASSERT_ALL_EQUAL( h(i,j), g(i,j) * 11 );
	}
	
	void testDivEquals() {
		*_g /= *_f;
		ASSERT_ALL_EQUAL( (*_g)(i,j), g(i,j) / f(i,j) );
	}

	void testDivEqualsDouble() {
		*_g /= 3;
		ASSERT_ALL_EQUAL( (*_g)(i,j), g(i,j) / 3 );
	}

	void testScalarDivScalar() {
		Scalar h = (*_g) / (*_f);
		ASSERT_ALL_EQUAL( h(i,j), g(i,j) / f(i,j) );
	}
	
	void testScalarDivDouble() {
		Scalar h = (*_g) / 11;
		ASSERT_ALL_EQUAL( h(i,j), g(i,j) / 11 );
	}
	
	void testDoubleDivScalar() {
		Scalar h = 11 / (*_g);
		ASSERT_ALL_EQUAL( h(i,j), 11 / g(i,j) );
	}
	
	void testUnaryMinus() {
		Scalar h = -(*_f);
		ASSERT_ALL_EQUAL( h(i,j), -f(i,j) );
	}
	
    void testIteratorCount() {
        Scalar::iterator i;
        int n=0;
        for (i = (*_f).begin(); i != (*_f).end(); ++i) {
            ++n;
        }
        TS_ASSERT_EQUALS(n, _nx * _ny);
    }

    void testIteratorAccess() {
        Scalar::iterator iter;
        int i=0;
        int j=0;
        for (iter = (*_f).begin(); iter != (*_f).end(); ++iter) {
            TS_ASSERT_DELTA( *iter, f(i,j), _delta );
            ++j;
            if (j >= _ny) {
                j -= _ny;
                ++i;
            }
        }
    }

    void testIteratorAssignment() {
        Scalar::iterator iter;
        int i=0;
        int j=0;
        for (iter = (*_f).begin(); iter != (*_f).end(); ++iter) {
            *iter = f(i,j) * 2 + 3;
            ++j;
            if (j >= _ny) {
                j -= _ny;
                ++i;
            }
        }
        ASSERT_ALL_EQUAL( (*_f)(i,j), f(i,j) * 2 + 3 );
    }

	void tearDown() {
		delete _f;
        delete _g;
        delete _grid;
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

#undef ASSERT_ALL_EQUAL