#include "Grid.h"
#include "Scalar.h"
#include "Flux.h"
#include "VectorOperations.h"
#include <gtest/gtest.h>

namespace {

class VectorOperationsTest : public testing::Test {
protected:
    VectorOperationsTest() : 
        _nx(3),
        _ny(5),
        _grid(_nx, _ny, 2, -1, -3),
        _f(_grid),
		_g(_grid),
		_p(_grid),
        _q(_grid)
    {
        // Initialize Scalar _f
		for (int i=0; i<_nx+1; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				_f(i,j) = f(i,j);
				_g(i,j) = g(i,j);				
			}
		}
		
		// Initialize Flux _q
		for (int i=0; i<_nx+1; ++i) {
			for (int j=0; j<_ny; ++j) {
				_q(X,i,j) = qx(i,j);
				_p(X,i,j) = px(i,j);				
			}
		}
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				_q(Y,i,j) = qy(i,j);
				_p(Y,i,j) = py(i,j);
			}
		}
    }
    
    // functions to test
	inline double f(int i, int j) {
		return 0.5 * i * i * _nx + 2 * j * (i + 1) + cos(j) * (_ny + 1);
	}

	inline double g(int i, int j) {
		return -5 * i * i + 3 * i * j + 2 * j * j;
	}
		
	inline double qx(int i, int j) {
		return 3 * i * i * _nx + 5 * (j +1 )* (1 + i + 0.7 * j) * cos(j);
	}
	
	inline double qy(int i, int j) {
		return  i * (i + j ) *3  + (_ny + cos(i)) * (j + 1 ) *2;
	}
	
	inline double px(int i, int j) {
		return 2*i + 3*j;
	}

	inline double py(int i, int j) {
		return i*j - 10;
	}
	
    // data
	int _nx;
	int _ny;
	Grid _grid;
	Scalar _f;
	Scalar _g;
	Flux _p;
	Flux _q;
};

// for Scalars
#define EXPECT_ALL_EQUAL(a,b)               \
	for (int i=0; i<_nx+1; ++i) {           \
		for (int j=0; j<_ny+1; ++j) {       \
			EXPECT_DOUBLE_EQ( (a), (b) );   \
		}                                   \
	}

// for Fluxes
#define EXPECT_ALL_X_EQUAL(a,b)             \
	for ( int i=0; i<_nx+1; ++i ) {         \
		for ( int j=0; j<_ny; ++j ) {       \
			EXPECT_DOUBLE_EQ( (a), (b) );   \
		}                                   \
	}

// for Fluxes
#define EXPECT_ALL_Y_EQUAL(a,b)             \
	for ( int i=0; i<_nx; ++i ) {           \
		for ( int j=0; j<_ny+1; ++j ) {     \
			EXPECT_DOUBLE_EQ( (a), (b) );   \
		}                                   \
	}

TEST_F(VectorOperationsTest, ScalarDotProductSymmetric) {
	EXPECT_DOUBLE_EQ(InnerProduct(_f,_g), InnerProduct(_f,_g) );        
}

TEST_F(VectorOperationsTest, ScalarDotProductDomainArea) {
	Scalar h( _grid );
	Scalar l( _grid );
	for (int i=0; i<_nx+1; ++i) {
		for (int j=0; j<_ny+1; ++j) {
			h(i,j) = 1;
			l(i,j) = 1;
		}
	}
	double area = _grid.getArea();
	EXPECT_DOUBLE_EQ(InnerProduct(h, l), area);
}

TEST_F(VectorOperationsTest, FluxDotProductSymmetric) {
	EXPECT_DOUBLE_EQ( InnerProduct(_q,_p), InnerProduct(_p, _q) );
}

TEST_F(VectorOperationsTest, FluxDotProductDomainArea) {
	Flux h( _grid );
	Flux l( _grid );
	for (int i=0; i<_nx+1; ++i) {
		for (int j=0; j<_ny; ++j) {
			h(X,i,j) = 1;
			l(X,i,j) = 0.2;
		}
	}
	for (int i=0; i<_nx; ++i) {
		for (int j=0; j<_ny+1; ++j) {
			h(Y,i,j) = 0.8;
			l(Y,i,j) = 1;
		}
	}
	double area = (_grid).getArea();
	EXPECT_DOUBLE_EQ(InnerProduct(h,l), area);
}

//	Tests on Scalar functions (that return Scalar):
TEST_F(VectorOperationsTest, ScalarCurlOfConstantEqualsZero) {		
	// Test 1: Curl of a non-zero constant Flux object
	// (w/ different Flux.x and Flux.y).
	Flux fluxc(_grid);
	double cx = 8;
	double cy = 5;
	fluxc = cx;
	for (int i=0; i<_nx; ++i) {
		for (int j=0; j<_ny+1; ++j) {
			fluxc(Y,i,j) = cy;
		}
	}
	
	Scalar sf = curl(fluxc);		
	EXPECT_ALL_EQUAL(sf(i,j), 0); 		
}

TEST_F(VectorOperationsTest, ScalarCurlOfFluxGivesZeroBCs) {		
	Scalar sf = curl(_q);		
	for (int i=0; i<_nx+1; ++i) {
		EXPECT_DOUBLE_EQ( sf(i,0), 0 );
		EXPECT_DOUBLE_EQ( sf(i,_ny), 0 );
	}
	for (int j=0; j<_ny+1; ++j) {       
		EXPECT_DOUBLE_EQ( sf(0,j), 0 );
		EXPECT_DOUBLE_EQ( sf(_nx,j), 0 );  
	}                         
}	

TEST_F(VectorOperationsTest, ScalarCurlOfSimplePolyFunctions) {
	Scalar sf = curl(_p); //p's x flux is 2i+3j; y flux is i*j-10.	
	for (int i = 1; i<_nx; ++i) {
		for (int j = 1; j<_ny; ++j) {
			EXPECT_DOUBLE_EQ( sf(i,j), -3+j ); 
		}	
	}	
	for (int i=0; i<_nx+1; ++i) {
		EXPECT_DOUBLE_EQ( sf(i,0), 0 );
		EXPECT_DOUBLE_EQ( sf(i,_ny), 0 );
	}
	for (int j=0; j<_ny+1; ++j) {       
		EXPECT_DOUBLE_EQ( sf(0,j), 0 );
		EXPECT_DOUBLE_EQ( sf(_nx,j), 0 );  
	}  		
}

//	Tests on Flux functions (that return Flux)
TEST_F(VectorOperationsTest, FluxCurlOfConstantEqualsZero) {
	Scalar scalarc(_grid);
	double c = 8;
	scalarc = c;
	Flux fluxf = curl(scalarc);	
	EXPECT_ALL_X_EQUAL(fluxf(X, i, j), 0);
	EXPECT_ALL_Y_EQUAL(fluxf(Y, i, j), 0);  	
}

TEST_F(VectorOperationsTest, FluxCurlOfSimplePolyFunction) {
	Flux fluxf = curl(_g);	// scalar g(i, j) = -5i^2+3*i*j+ 2j^2;
	EXPECT_ALL_X_EQUAL(fluxf(X, i, j), 3*i+4*j+2);
	EXPECT_ALL_Y_EQUAL(fluxf(Y, i, j), 10*i-3*j+5);  	
}


#undef EXPECT_ALL_EQUAL
#undef EXPECT_ALL_X_EQUAL
#undef EXPECT_ALL_Y_EQUAL

} // namespace