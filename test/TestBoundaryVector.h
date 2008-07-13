#include <cxxtest/TestSuite.h>
#include "BoundaryVector.h"

#define ASSERT_ALL_EQUAL(a,b) \
	for ( int i=0; i<_n; ++i ) { \
		TS_ASSERT_DELTA( (a), (b), _delta ); \
	}
	
class TestBoundaryVector : public CxxTest::TestSuite {
public:
	void setUp() {
		_n = 10;
		_f = new BoundaryVector( _n );
		_g = new BoundaryVector( _n );
		for (int i=0; i<_n; ++i) {
			_f->x(i) = fx(i);
			_f->y(i) = fy(i);
			_g->x(i) = gx(i);
			_g->y(i) = gy(i);
		}
	}

	void testAssignment() {
		int i = _n/2;
		_f->x(i) = 73;
		TS_ASSERT_DELTA(_f->x(i), 73, _delta);
		_f->y(i) = 23;
		TS_ASSERT_DELTA(_f->y(i), 23, _delta);
	}

	void testCopyConstructor() {
		BoundaryVector h = *_f;
		ASSERT_ALL_EQUAL( h.x(i), fx(i) );
		ASSERT_ALL_EQUAL( h.y(i), fy(i) );
	}

	void testCopyFromDouble() {
		double a = 7;
		*_f = a;
		ASSERT_ALL_EQUAL( _f->x(i), a );
		ASSERT_ALL_EQUAL( _f->y(i), a );
	}

	void testUnaryMinus() {
		BoundaryVector h = -(*_f);
		ASSERT_ALL_EQUAL( h.x(i), -fx(i) );
		ASSERT_ALL_EQUAL( h.y(i), -fy(i) );
	}

	void testPlusEqualsBoundaryVector() {
		*_g += *_f;
		ASSERT_ALL_EQUAL( _g->x(i), fx(i) + gx(i) );
		ASSERT_ALL_EQUAL( _g->y(i), fy(i) + gy(i) );
	}
	
	void testMinusEqualsBoundaryVector() {
		*_g -= *_f;
		ASSERT_ALL_EQUAL( _g->x(i), gx(i) - fx(i) );
		ASSERT_ALL_EQUAL( _g->y(i), gy(i) - fy(i) );
	}
	
	void testBoundaryVectorPlusBoundaryVector() {
		BoundaryVector h = (*_f) + (*_g);
		ASSERT_ALL_EQUAL( h.x(i), fx(i) + gx(i) );
		ASSERT_ALL_EQUAL( h.y(i), fy(i) + gy(i) );
	}
	
	void testBoundaryVectorMinusBoundaryVector() {
		BoundaryVector h = (*_f) - (*_g);
		ASSERT_ALL_EQUAL( h.x(i), fx(i) - gx(i) );
		ASSERT_ALL_EQUAL( h.y(i), fy(i) - gy(i) );
	}
	
	void testTimesEqualsDouble() {
		(*_f) *= 3;
		ASSERT_ALL_EQUAL( _f->x(i), 3 * fx(i) );
		ASSERT_ALL_EQUAL( _f->y(i), 3 * fy(i) );
	}

	void testDivEqualsDouble() {
		(*_f) /= 3.;
		ASSERT_ALL_EQUAL( _f->x(i), fx(i) / 3. );
		ASSERT_ALL_EQUAL( _f->y(i), fy(i) / 3. );
	}
	
	void testTimesDouble() {
		double a = 7;
		BoundaryVector h = (*_f) * a;
		ASSERT_ALL_EQUAL( h.x(i), fx(i) * a );
		ASSERT_ALL_EQUAL( h.y(i), fy(i) * a );
	}

	void testDoubleTimes() {
		double a = 7;
		BoundaryVector h = a * (*_f);
		ASSERT_ALL_EQUAL( h.x(i), fx(i) * a );
		ASSERT_ALL_EQUAL( h.y(i), fy(i) * a );
	}

	void testDivDouble() {
		double a = 7;
		BoundaryVector h = (*_f) / a;
		ASSERT_ALL_EQUAL( h.x(i), fx(i) / a );
		ASSERT_ALL_EQUAL( h.y(i), fy(i) / a );
	}

    void testIteratorCount() {
        BoundaryVector::iterator i;
        int n;

        // Loop over X and Y dimensions
        for (int dim=BoundaryVector::X; dim <= BoundaryVector::Y; ++dim) {
            // Count number of elements
            n = 0;
            for (i = (*_f).begin(dim); i != (*_f).end(dim); ++i) {
                ++n;
            }
            TS_ASSERT_EQUALS(n, _n);
        }
    }

    void testIteratorAccess() {
        BoundaryVector::iterator iter;
        int i=0;

        // X dimension
        for (iter=_f->begin(BoundaryVector::X); iter != _f->end(BoundaryVector::X); ++iter) {
            TS_ASSERT_DELTA( *iter, fx(i), _delta );
            ++i;
        }
        
        // Y dimension
        i = 0;
        for (iter=_f->begin(BoundaryVector::Y); iter != _f->end(BoundaryVector::Y); ++iter) {
            TS_ASSERT_DELTA( *iter, fy(i), _delta );
            ++i;
        }
    }

    void testIteratorAssignment() {
        BoundaryVector::iterator iter;
        int i=0;

        // X dimension
        for (iter=_f->begin(BoundaryVector::X); iter != _f->end(BoundaryVector::X); ++iter) {
            *iter = 2 * fx(i) + 3;
            ++i;
        }
        ASSERT_ALL_EQUAL( _f->x(i), 2 * fx(i) + 3 );
        ASSERT_ALL_EQUAL( _f->y(i), fy(i) ); // y-dir unchanged

        // Y dimension
        i = 0;
        for (iter=_g->begin(BoundaryVector::Y); iter != _g->end(BoundaryVector::Y); ++iter) {
            *iter = 3 * gy(i) + 2;
            ++i;
        }
        ASSERT_ALL_EQUAL( _g->x(i), gx(i) );  // x-dir unchanged
        ASSERT_ALL_EQUAL( _g->y(i), 3 * gy(i) + 2 );
    }

	void tearDown() {
		delete _f;
        delete _g;
	}

	double fx(int i) {
		return 2*i;
	}
	
	double fy(int i) {
		return 3*i;
	}

	double gx(int i) {
		return _n-i;
	}

	double gy(int i) {
		return 2*_n-i;
	}

private:
	int _n;
	BoundaryVector* _f;
	BoundaryVector* _g;
	static double _delta;
};

double TestBoundaryVector::_delta = 1e-10;

#undef ASSERT_ALL_EQUAL