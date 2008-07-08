#include <cxxtest/TestSuite.h>
#include "Grid.h"
#include "Flux.h"

#define ASSERT_ALL_X_EQUAL(a,b) \
	for ( int i=0; i<_nx+1; ++i ) { \
		for ( int j=0; j<_ny; ++j ) { \
			TS_ASSERT_DELTA( (a), (b), _delta ); \
		} \
	}
#define ASSERT_ALL_Y_EQUAL(a,b) \
	for ( int i=0; i<_nx; ++i ) { \
		for ( int j=0; j<_ny+1; ++j ) { \
			TS_ASSERT_DELTA( (a), (b), _delta ); \
		} \
	}
	
class TestFlux : public CxxTest::TestSuite {
public:
	void setUp() {
		// define a Grid
		_nx = 3;
		_ny = 5;
		double length = 2;
		double xOffset = -1;
		double yOffset = -3;
		_grid = new Grid( _nx, _ny, length, xOffset, yOffset );

		_f = new Flux( *_grid );
		_g = new Flux( *_grid );
		for (int i=0; i<_nx+1; ++i) {
			for (int j=0; j<_ny; ++j) {
				(*_f).x(i,j) = fx(i,j);
				(*_g).x(i,j) = gx(i,j);
			}
		}
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				(*_f).y(i,j) = fy(i,j);
				(*_g).y(i,j) = gy(i,j);
			}
		}
	}

	void testAssignment() {
		_f->x(1,3) = 73;
		TS_ASSERT_DELTA(_f->x(1,3), 73, _delta);
	}

	void testCopyConstructor() {
		Flux h = *_f;
		ASSERT_ALL_X_EQUAL( h.x(i,j), fx(i,j) );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), fy(i,j) );
	}

	void testCopyFromDouble() {
		double a = 7;
		*_f = a;
		ASSERT_ALL_X_EQUAL( _f->x(i,j), a );
		ASSERT_ALL_Y_EQUAL( _f->y(i,j), a );
	}

	void testUnaryMinus() {
		Flux h = -(*_f);
		ASSERT_ALL_X_EQUAL( h.x(i,j), -fx(i,j) );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), -fy(i,j) );
	}

	void testPlusEqualsFlux() {
		*_g += *_f;
		ASSERT_ALL_X_EQUAL( _g->x(i,j), fx(i,j) + gx(i,j) );
		ASSERT_ALL_Y_EQUAL( _g->y(i,j), fy(i,j) + gy(i,j) );
	}
	
	void testPlusEqualsDouble() {
		*_g += 7;
		ASSERT_ALL_X_EQUAL( _g->x(i,j), gx(i,j) + 7 );
		ASSERT_ALL_Y_EQUAL( _g->y(i,j), gy(i,j) + 7 );
	}
	
	void testMinusEqualsFlux() {
		*_g -= *_f;
		ASSERT_ALL_X_EQUAL( _g->x(i,j), gx(i,j) - fx(i,j) );
		ASSERT_ALL_Y_EQUAL( _g->y(i,j), gy(i,j) - fy(i,j) );
	}
	
	void testMinusEqualsDouble() {
		*_g -= 7;
		ASSERT_ALL_X_EQUAL( _g->x(i,j), gx(i,j) - 7 );
		ASSERT_ALL_Y_EQUAL( _g->y(i,j), gy(i,j) - 7 );
	}
	
	void testFluxPlusFlux() {
		Flux h = (*_f) + (*_g);
		ASSERT_ALL_X_EQUAL( h.x(i,j), fx(i,j) + gx(i,j) );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), fy(i,j) + gy(i,j) );
	}
	
	void testFluxPlusDouble() {
		Flux h = (*_f) + 7;
		ASSERT_ALL_X_EQUAL( h.x(i,j), fx(i,j) + 7 );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), fy(i,j) + 7 );		
	}

	void testDoublePlusFlux() {
		Flux h = 7 + (*_f);
		ASSERT_ALL_X_EQUAL( h.x(i,j), fx(i,j) + 7 );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), fy(i,j) + 7 );		
	}

	void testFluxMinusFlux() {
		Flux h = (*_f) - (*_g);
		ASSERT_ALL_X_EQUAL( h.x(i,j), fx(i,j) - gx(i,j) );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), fy(i,j) - gy(i,j) );
	}
	
	void testFluxMinusDouble() {
		Flux h = (*_f) - 7;
		ASSERT_ALL_X_EQUAL( h.x(i,j), fx(i,j) - 7 );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), fy(i,j) - 7 );		
	}

	void testDoubleMinusFlux() {
		Flux h = 7 - (*_f);
		ASSERT_ALL_X_EQUAL( h.x(i,j), 7 - fx(i,j) );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), 7 - fy(i,j) );		
	}
	
	void testTimesEqualsDouble() {
		(*_f) *= 3;
		ASSERT_ALL_X_EQUAL( _f->x(i,j), 3 * fx(i,j) );
		ASSERT_ALL_Y_EQUAL( _f->y(i,j), 3 * fy(i,j) );
	}

	void testDivEqualsDouble() {
		(*_f) /= 3.;
		ASSERT_ALL_X_EQUAL( _f->x(i,j), fx(i,j) / 3. );
		ASSERT_ALL_Y_EQUAL( _f->y(i,j), fy(i,j) / 3. );
	}

	void testFluxTimesDouble() {
		double a = 5;
		Flux h = (*_f) * a;
		ASSERT_ALL_X_EQUAL( h.x(i,j), fx(i,j) * a );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), fy(i,j) * a );		
	}

	void testDoubleTimesFlux() {
		double a = 5;
		Flux h = a * (*_f);
		ASSERT_ALL_X_EQUAL( h.x(i,j), fx(i,j) * a );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), fy(i,j) * a );		
	}

	void testFluxDivDouble() {
		double a = 5;
		Flux h = (*_f) / a;
		ASSERT_ALL_X_EQUAL( h.x(i,j), fx(i,j) / a );
		ASSERT_ALL_Y_EQUAL( h.y(i,j), fy(i,j) / a );		
	}

	void tearDown() {
		delete _f;
	}

	double fx(int i, int j) {
		return i * _nx + j;
	}
	
	double fy(int i, int j) {
		return j * _ny - i;
	}

	double gx(int i, int j) {
		return i + 3*j;
	}

	double gy(int i, int j) {
		return i*j - 10;
	}

private:
	int _nx;
	int _ny;
	Grid* _grid;
	Flux* _f;
	Flux* _g;
	static double _delta;
};

double TestFlux::_delta = 1e-10;

#undef ASSERT_ALL_X_EQUAL
#undef ASSERT_ALL_Y_EQUAL
