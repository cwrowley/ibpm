#include <cxxtest/TestSuite.h>
#include "Grid.h"
#include "Scalar.h"
#include "Flux.h"
#include "VectorOperations.h"

#define ASSERT_ALL_EQUAL(a,b)                        \
	for (int i=0; i<_nx+1; ++i) {                      \
		for (int j=0; j<_ny+1; ++j) {                  \
			TS_ASSERT_DELTA((a), (b), _delta);       \
		}                                            \
	}
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


class TestVectorOperations : public CxxTest::TestSuite {
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
		for (int i=0; i<_nx+1; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				(*_f)(i,j) = f(i,j);				
			}
		}
		
		_q= new Flux( *_grid );
		for (int i=0; i<_nx+1; ++i) {
			for (int j=0; j<_ny; ++j) {
				(*_q).x(i,j) = qx(i,j);				
			}
		}
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				(*_q).y(i,j) = qy(i,j);
			}
		}
	}	

//	Tests on Scalar functions (that return Scalar):
   void testScalarCurlOfConstantEqualsZero() {		
		// Test 1: Curl of a non-zero constant Flux object (w/ different Flux.x and Flux.y).
		Flux fluxc(*_grid);
		double cx = 8;
		double cy = 5;
		fluxc = cx;
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				fluxc.y(i,j) = cy;
			}
		}
		
		Scalar sf = curl(fluxc);		
		ASSERT_ALL_EQUAL(sf(i,j), 0); 		
    }

	void testScalarCurlValue() {
		// Test 2: Curl of a non-constant Flux object.			
		Scalar h(*_grid);
		for (int i = 1; i < _nx; ++i) {
			for (int j=1; j<_ny; ++j) {
				h(i,j) = qy(i,j) - qy(i-1,j) - qx(i,j) + qx(i,j-1);
			}			
		}
		for (int j = 0; j < _ny+1; ++j) {
			h(0,j) = 0;
			h(_nx,j) = 0;			
		}
		for (int i = 0; i < _nx+1; ++i) {
			h(i,0) = 0;
			h(i,_ny) = 0;			
		}

		Scalar sf = curl(*_q);		
		ASSERT_ALL_EQUAL(sf(i,j), h(i,j)); 
	}
	
//	void testCurlGradientEqualsZero() {
//		// Assume zero b.cs at the bottom and left "ghost edges".
//		Flux fluxgrad(*_grid);
//		fluxgrad.gradient(*_f);
//		_f->curl(fluxgrad);
//		ASSERT_ALL_EQUAL( (*_f)(i,j), 0 );
//    }
//    
//	void testDivCurlEqualsZero() {
//		// Divergence of Curl of a non-constant Scalar object.
//		Flux fluxcurl(*_grid);
//		fluxcurl.curl(*_f);		
//		(*_f).divergence(fluxcurl);		
//		ASSERT_ALL_EQUAL( (*_f)(i,j), 0 );		
//    }
//		
//	void testDivConstantEqualsZero() {
//		// Divergence of a non-zero constant flux (w/ different Flux.x and Flux.y)
//		Flux fluxc(*_grid);
//		double cx = 8;
//		double cy = 7;
//		fluxc = cx;
//		for (int i=0; i<_nx; ++i) {
//			for (int j=0; j<_ny+1; ++j) {
//				fluxc.y(i,j) = cy;
//			}
//		}
//		(*_f).divergence(fluxc);
//		ASSERT_ALL_EQUAL((*_f)(i,j), 0);		
//    }
//    
//	void testDivNonConstantFlux() {
//		// Divergence of a non-constant flux.
//		Flux fluxf(*_grid);
//		for (int i=0; i<_nx+1; ++i) {
//			for (int j=0; j<_ny; ++j) {
//				fluxf.x(i,j) = fx(i,j);				
//			}
//		}
//		for (int i=0; i<_nx; ++i) {
//			for (int j=0; j<_ny+1; ++j) {
//				fluxf.y(i,j) = fy(i,j);				
//			}
//		}		
//
//		(*_f).divergence(fluxf);
//		ASSERT_ALL_EQUAL( (*_f)(i,j), fx(i+1,j) - fx(i,j) + fy(i,j+1) -fy(i,j) );	
//	}	
//

//	Tests on Flux functions (that return Flux)
	void testFluxCurlOfConstantEqualsZero() {
		Scalar scalarc(*_grid);
		double c = 8;
		scalarc = c;
		Flux fluxf = curl(scalarc);	
		ASSERT_ALL_X_EQUAL(fluxf.x(i, j), 0);
		ASSERT_ALL_Y_EQUAL(fluxf.y(i, j), 0);  	
    }

	void testCurlValue() {
		// Curl of a non-constant Scalar object.				
		Flux h(*_grid);
		// X direction.
		for (int i = 0; i < _nx+1; ++i) {
			for ( int j = 0; j < _ny; ++j) {
				h.x(i,j) = f(i,j+1) - f(i,j);
			}					
		}		
			
		// Y direction.
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				h.y(i,j) = f(i,j) - f(i+1,j);
			}									
		}
		
		Flux ff=curl(*_f);		
		ASSERT_ALL_X_EQUAL(ff.x(i,j), h.x(i,j)); 
		ASSERT_ALL_Y_EQUAL(ff.y(i,j), h.y(i,j));		
	}
//	
//	void testGradientOfConstantEqualsZero() {
//		// Assume zero b.cs at the ghost cells (-1,:), (nx, :), (:,-1), (:, ny).
//		// Gradient of a non-zero constant Scalar object.
//		// Note: A test of curl(gradient(Scalar object) is given in
//		//       TestScalar.h.)		
//		Scalar scalarc(*_grid);
//		double c = 6;
//		scalarc = c;
//		(*_f).gradient(scalarc);
//		// Check inner grids.
//		for ( int i=1; i<_nx; ++i ) { 
//			for ( int j=0; j<_ny; ++j ) { 
//				TS_ASSERT_DELTA( (*_f).x(i,j), 0, _delta ); 
//			} 
//		}
//		for ( int i=0; i<_nx; ++i ) { 
//			for ( int j=1; j<_ny; ++j ) { 
//				TS_ASSERT_DELTA( (*_f).y(i,j), 0, _delta ); 
//			} 
//		}
//		// Check boundary grids.
//		// CWR: ???  Why should boundary values equal the constant?
//		for ( int j=0; j<_ny; ++j ) { 
//			TS_ASSERT_DELTA( (*_f).x(0,j), c, _delta );
//			TS_ASSERT_DELTA( (*_f).x(_nx,j), -c, _delta );  
//		}
//		for ( int i=0; i<_nx; ++i ) { 
//			TS_ASSERT_DELTA( (*_f).y(i,0), c, _delta );
//			TS_ASSERT_DELTA( (*_f).y(i,_ny), -c, _delta );
//		}
//    }
//
//	void testGradientValue() {
//		// Test 2: Gradient of a non-constant Scalar object.
//		Scalar scalarf(*_grid);
//		for (int i=0; i<_nx; ++i) {
//			for (int j=0; j<_ny; ++j) {
//				scalarf(i,j) = f(i,j);				
//			}
//		}
//		
//		Flux h( *_grid );
//		// X direction.		
//		for ( int j = 0; j < _ny; ++j) {
//				h.x(0,j) = f(0,j);
//				h.x(_nx,j) = -f(_nx-1,j);
//			}
//		for (int i = 1; i < _nx; ++i) {
//			for ( int j = 0; j < _ny; ++j) {
//				h.x(i,j) = f(i,j) - f(i-1,j);
//			}			
//		}
//		// Y direction.
//		for ( int i = 0; i < _nx; ++i) {
//			h.y(i,0) = f(i,0);
//			h.y(i,_ny) = -f(i,_ny-1);
//		}
//		for (int i=0; i<_nx; ++i) {
//			for (int j=1; j<_ny; ++j) {
//				h.y(i,j) = f(i,j) - f(i,j-1);
//			}			
//		}
//		(*_f).gradient(scalarf); 		
//		ASSERT_ALL_X_EQUAL((*_f).x(i,j), h.x(i,j)); 
//		ASSERT_ALL_Y_EQUAL((*_f).y(i,j), h.y(i,j));
//	}
	
	
	void tearDown() {
		delete _f;
        delete _q;
        delete _grid;
	}

	// functions to test
	double f(int i, int j) {
		return 0.5 * i * i * _nx + 2 * j * (i + 1) + cos(j) * (_ny + 1);
	}

	double g(int i, int j) {
		return 3*i + j;
	}
	
	double qx(int i, int j) {
		return 3 * i * i * _nx + 5 * (j +1 )* (1 + i + 0.7 * j) * cos(j);
	}
	
	double qy(int i, int j) {
		return  i * (i + j ) *3  + (_ny + cos(i)) * (j + 1 ) *2;
	}

private:
	int _nx;
	int _ny;
	Scalar* _f;
	Flux* _q;
	Grid* _grid;
	static double _delta;
};

double TestVectorOperations::_delta = 1e-10;

#undef ASSERT_ALL_EQUAL
#undef ASSERT_ALL_X_EQUAL
#undef ASSERT_ALL_Y_EQUAL