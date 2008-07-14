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
	
	void testCurlOfConstantEqualsZero() {
		// Assume zero b.cs at the ghost cells (nx, :), (:, ny).
		Scalar scalarc(*_grid);
		double c = 8;
		scalarc = c;
		(*_f).curl(scalarc);		
		// X direction.
		for ( int i=0; i<_nx; ++i ) { 
			for ( int j=0; j<_ny-1; ++j ) { 
				TS_ASSERT_DELTA( (*_f).x(i,j), 0, _delta );				 
			}
			TS_ASSERT_DELTA( (*_f).x(i,_ny-1), -c, _delta ); 
		}
		for ( int j=0; j<_ny; ++j ) { 
				TS_ASSERT_DELTA( (*_f).x(_nx,j), 0, _delta );				 
		}
		// Y direction.
		for ( int j=0; j<_ny; ++j ) { 
			for ( int i=0; i<_nx-1; ++i ) {			
				TS_ASSERT_DELTA( (*_f).y(i,j), 0, _delta ); 
			}
			TS_ASSERT_DELTA( (*_f).y(_nx-1,j), c, _delta ); 
		}
		for ( int i=0; i<_nx; ++i ) {			
				TS_ASSERT_DELTA( (*_f).y(i,_ny), 0, _delta ); 
		}
    }

	void testCurlValue() {
		// Curl of a non-constant Scalar object. 		
		Scalar scalarf(*_grid);		
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny; ++j) {
				scalarf(i,j) = f(i,j);				
			}
		}
		
		Flux h(*_grid);
		// X direction.
		for (int i = 0; i < _nx; ++i) {
			for ( int j = 0; j < _ny-1; ++j) {
				h.x(i,j) = f(i,j+1) - f(i,j);
			}
			h.x(i,_ny-1) = - f (i,_ny-1);			
		}
		for (int j = 0; j < _ny; ++j) {
			h.x(_nx, j) = 0;
		}		
			
		// Y direction.
		for (int j=0; j<_ny; ++j) {
			for (int i=0; i<_nx-1; ++i) {
				h.y(i,j) = f(i,j) - f(i+1,j);
			}
			h.y(_nx-1,j) = f(_nx-1,j);						
		}
		for (int i = 0; i < _nx; ++i) {
			h.y(i, _ny) = 0;
		}
		(*_f).curl(scalarf);		
		ASSERT_ALL_X_EQUAL((*_f).x(i,j), h.x(i,j)); 
		ASSERT_ALL_Y_EQUAL((*_f).y(i,j), h.y(i,j));		
	}
	
	void testGradientOfConstantEqualsZero() {
		// Assume zero b.cs at the ghost cells (-1,:), (nx, :), (:,-1), (:, ny).
		// Gradient of a non-zero constant Scalar object.
		// Note: A test of curl(gradient(Scalar object) is given in
		//       TestScalar.h.)		
		Scalar scalarc(*_grid);
		double c = 6;
		scalarc = c;
		(*_f).gradient(scalarc);
		// Check inner grids.
		for ( int i=1; i<_nx; ++i ) { 
			for ( int j=0; j<_ny; ++j ) { 
				TS_ASSERT_DELTA( (*_f).x(i,j), 0, _delta ); 
			} 
		}
		for ( int i=0; i<_nx; ++i ) { 
			for ( int j=1; j<_ny; ++j ) { 
				TS_ASSERT_DELTA( (*_f).y(i,j), 0, _delta ); 
			} 
		}
		// Check boundary grids.
		// CWR: ???  Why should boundary values equal the constant?
		for ( int j=0; j<_ny; ++j ) { 
			TS_ASSERT_DELTA( (*_f).x(0,j), c, _delta );
			TS_ASSERT_DELTA( (*_f).x(_nx,j), -c, _delta );  
		}
		for ( int i=0; i<_nx; ++i ) { 
			TS_ASSERT_DELTA( (*_f).y(i,0), c, _delta );
			TS_ASSERT_DELTA( (*_f).y(i,_ny), -c, _delta );
		}
    }

	void testGradientValue() {
		// Test 2: Gradient of a non-constant Scalar object.
		Scalar scalarf(*_grid);
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny; ++j) {
				scalarf(i,j) = f(i,j);				
			}
		}
		
		Flux h( *_grid );
		// X direction.		
		for ( int j = 0; j < _ny; ++j) {
				h.x(0,j) = f(0,j);
				h.x(_nx,j) = -f(_nx-1,j);
			}
		for (int i = 1; i < _nx; ++i) {
			for ( int j = 0; j < _ny; ++j) {
				h.x(i,j) = f(i,j) - f(i-1,j);
			}			
		}
		// Y direction.
		for ( int i = 0; i < _nx; ++i) {
			h.y(i,0) = f(i,0);
			h.y(i,_ny) = -f(i,_ny-1);
		}
		for (int i=0; i<_nx; ++i) {
			for (int j=1; j<_ny; ++j) {
				h.y(i,j) = f(i,j) - f(i,j-1);
			}			
		}
		(*_f).gradient(scalarf); 		
		ASSERT_ALL_X_EQUAL((*_f).x(i,j), h.x(i,j)); 
		ASSERT_ALL_Y_EQUAL((*_f).y(i,j), h.y(i,j));
	}
	
    void testDotProductSymmetric() {
		TS_ASSERT_DELTA( _f->dot(*_g), _g->dot(*_f), _delta );
    }

	void testDotProductValue() {
		double dp = 0;
		for (int i = 0; i < _nx+1; ++i) {
			for ( int j = 0; j < _ny; ++j) {
				dp += fx(i, j) * gx(i, j);
			}			
		}
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				dp += fy(i, j) * gy(i, j);
			}			
		}
		TS_ASSERT_DELTA( _f->dot(*_g), dp, _delta ); 
	}

    void testIteratorCount() {
        Flux::iterator i;
        int n;

        // X direction
        // Number of fluxes in X-dir is (nx+1, ny)
        n = 0;
        for (i = _f->begin(Flux::X); i != _f->end(Flux::X); ++i) {
            ++n;
        }
        TS_ASSERT_EQUALS(n, (_nx + 1) * _ny );

        // Y direction
        // Number of fluxes in Y-dir is (nx, ny+1)
        n = 0;
        for (i = _f->begin(Flux::Y); i != _f->end(Flux::Y); ++i) {
            ++n;
        }
        TS_ASSERT_EQUALS(n, _nx * (_ny + 1) );
    }

    void testIteratorAccess() {
        Flux::iterator iter;
        int i=0;
        int j=0;
    
        // X dimension
        for (iter=_f->begin(Flux::X); iter != _f->end(Flux::X); ++iter) {
            TS_ASSERT_DELTA( *iter, fx(i,j), _delta );
            ++j;
            if (j >= _ny) {
                j -= _ny;
                ++i;
            }
        }
        
        // Y dimension
        // Caution: number of fluxes in Y-dir is ny+1 !
        i = 0;
        j = 0;
        for (iter=_f->begin(Flux::Y); iter != _f->end(Flux::Y); ++iter) {
            TS_ASSERT_DELTA( *iter, fy(i,j), _delta );
            ++j;
            if (j >= (_ny+1)) {
                j -= (_ny+1);
                ++i;
            }
        }
    }
    
    void testIteratorAssignment() {
        Flux::iterator iter;
        int i=0;
        int j=0;
    
        // X dimension
        for (iter=_f->begin(Flux::X); iter != _f->end(Flux::X); ++iter) {
            *iter = 2 * fx(i,j) + 3;
            ++j;
            if (j >= _ny) {
                j -= _ny;
                ++i;
            }
        }
        ASSERT_ALL_X_EQUAL( _f->x(i,j), 2 * fx(i,j) + 3 );
        ASSERT_ALL_Y_EQUAL( _f->y(i,j), fy(i,j) );

        // Y dimension
        // Caution: number of fluxes in Y-dir is ny+1 !
        i = 0;
        j = 0;
        for (iter=_g->begin(Flux::Y); iter != _g->end(Flux::Y); ++iter) {
            *iter = 3 * gy(i,j) + 2;
            ++j;
            if (j >= (_ny+1)) {
                j -= (_ny+1);
                ++i;
            }
        }
        ASSERT_ALL_X_EQUAL( _g->x(i,j), gx(i,j) );
        ASSERT_ALL_Y_EQUAL( _g->y(i,j), 3 * gy(i,j) + 2 );

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
	
	double f(int i, int j) {
		return 0.5 * i * i * _nx + 2 * j * (i + 1) + cos(j) * (_ny + 1);
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
