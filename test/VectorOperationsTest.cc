#include "Grid.h"
#include "Scalar.h"
#include "Flux.h"
#include "BoundaryVector.h"
#include "VectorOperations.h"
#include "SingleWavenumber.h"
#include <gtest/gtest.h>

using namespace ibpm;

namespace {

class VectorOperationsTest : public testing::Test {
protected:
    VectorOperationsTest() : 
        _nx(4),
        _ny(8),
        _grid(_nx, _ny, 2, -1, -3),
        _f(_grid),
		_g(_grid),
		_x(_grid),
		_y(_grid),
		_p(_grid),
        _q(_grid)
    {
        // Initialize Scalars _f, _g, _x, _y
		for (int i=0; i<_nx+1; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				_f(i,j) = f(i,j);
				_g(i,j) = g(i,j);		
                _x(i,j) = _grid.getXEdge(i);
                _y(i,j) = _grid.getYEdge(j);
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
    // TODO: Pick better functions, with analytic expressions (in x and y)
    // (Single wavenumber for the scalar fields?)
	inline double f(int i, int j) {
		return 0.5 * i * i * _nx + 2 * j * (i + 1) + cos(j) * (_ny + 1);
	}

	inline double g(int i, int j) {
		return -5 * i * i + 3 * i * j + 2 * j * j;
	}
		
	inline double qx(int i, int j) {
		return 3 * i * _nx + 5 * (j +1 )* i;
	}
	
	inline double qy(int i, int j) {
		return 4 * i *  j   + _ny + cos(i);
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
    Scalar _x;
    Scalar _y;
	Flux _p;
	Flux _q;
};

// for Scalars
#define EXPECT_ALL_EQ(a,b)                  \
	for (int i=0; i<_nx+1; ++i) {           \
		for (int j=0; j<_ny+1; ++j) {       \
			EXPECT_DOUBLE_EQ( (a), (b) );   \
		}                                   \
	}

// for Scalars
#define EXPECT_ALL_INTERIOR_EQ(a,b)         \
	for (int i=1; i<_nx; ++i) {             \
		for (int j=1; j<_ny; ++j) {         \
			EXPECT_DOUBLE_EQ( (a), (b) );   \
		}                                   \
	}

// for Fluxes
#define EXPECT_ALL_X_EQ(a,b)                \
	for ( int i=0; i<_nx+1; ++i ) {         \
		for ( int j=0; j<_ny; ++j ) {       \
			EXPECT_DOUBLE_EQ( (a), (b) );   \
		}                                   \
	}

// for Fluxes
#define EXPECT_ALL_Y_EQ(a,b)                \
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
	double area = _grid.getArea();
    double dx = _grid.getDx();
	EXPECT_DOUBLE_EQ(InnerProduct(h,l), area  / (dx * dx) );
}

// =============================
// = Flux to velocity and back =
// =============================

// < q, X^* u > = < X q, u > for all q and u
//    X is the map from from Flux (at edges) to Scalar x-velocity (at nodes)
//    X^* is its adjoint, a map from Scalar x-velocity to Flux
TEST_F(VectorOperationsTest, FluxToXVelocity) {
    Flux q1(_grid);
    Flux q2(_grid);
    Scalar u1(_grid);
    Scalar u2(_grid);
    
    // TODO: Eventually, loop over several different values of q1 and u2
    q1 = 2.;
    u2 = 5.;
    q2 = 0;
    XVelocityToFlux( u2, q2 );
    FluxToXVelocity( q1, u1 );
    double IPFlux = InnerProduct( q1, q2 );
    double IPScalar = InnerProduct( u1, u2 );
    EXPECT_DOUBLE_EQ( IPFlux, IPScalar );
}

// < q, Y^* v > = < Y q, v > for all q and v
//    Y is the map from from Flux (at edges) to Scalar y-velocity (at nodes)
//    Y^* is its adjoint, a map from Scalar y-velocity to Flux
TEST_F(VectorOperationsTest, FluxToYVelocity) {
    Flux q1(_grid);
    Flux q2(_grid);
    Scalar v1(_grid);
    Scalar v2(_grid);
    
    // TODO: Eventually, loop over several different values of q1 and v2
    q1 = 3.;
    v2 = 7.;
    q2 = 0;
    YVelocityToFlux( v2, q2 );
    FluxToYVelocity( q1, v1 );
    double IPFlux = InnerProduct( q1, q2 );
    double IPScalar = InnerProduct( v1, v2 );
    EXPECT_DOUBLE_EQ( IPFlux, IPScalar );
}

TEST_F(VectorOperationsTest, FluxToVelocity) {
    Flux q1(_grid);
    Flux q2(_grid);
    Scalar u1(_grid);
    Scalar u2(_grid);
    Scalar v1(_grid);
    Scalar v2(_grid);
    
    // TODO: Eventually, loop over several different values of q1 and (u2, v2)
    q1 = 3.;
    u2 = 5.;
    v2 = 7.;
    q2 = 0;
    VelocityToFlux( u2, v2, q2 );
    FluxToVelocity( q1, u1, v1 );
    double IPFlux = InnerProduct( q1, q2 );
    double IPScalar = InnerProduct( u1, u2 ) + InnerProduct( v1, v2 );
    EXPECT_DOUBLE_EQ( IPFlux, IPScalar );    
}

// =================
// = Cross product =
// =================

// Test that q x q = 0
TEST_F(VectorOperationsTest, CrossProductWithSelf) {
    Flux q(_grid);
    
    // TODO: Eventually, test several different values of q
    Scalar a = CrossProduct( q, q );
    EXPECT_ALL_EQ( 0, a(i,j) );
}

// Test that (u,0) x (0,v) = u v
TEST_F(VectorOperationsTest, CrossProductOfOrthogonalVectors) {
    Scalar u(_grid);
    Scalar v(_grid);
    Flux q1(_grid);
    Flux q2(_grid);

    // TODO: Everntually, test several different values of u, v
    u = 3;
    v = 5;
    q1 = 0;
    q2 = 0;
    XVelocityToFlux( u, q1 );
    YVelocityToFlux( v, q2 );
    Scalar cross = CrossProduct( q1, q2 );
    EXPECT_ALL_EQ( u(i,j)*v(i,j), cross(i,j) );
}

// Test that q1 x q2 = -q2 x q1
TEST_F(VectorOperationsTest, CrossProductAntiSymmetric) {
    Flux q1(_grid);
    Flux q2(_grid);
    
    // TODO: Eventually, test several different values of q1, q2
    q1 = 5;
    q2 = 3;
    Scalar q1q2 = CrossProduct( q1, q2 );
    Scalar q2q1 = CrossProduct( q2, q1 );
    EXPECT_ALL_EQ( q1q2(i,j), -q2q1(i,j) );
}

// Test that < a, q1 x q2 > = < q1, q2 x a >
// for all Scalars a, Fluxes q1, q2
TEST_F(VectorOperationsTest, CrossProductTripleProduct) {
    Flux q1(_grid);
    Flux q2(_grid);
    Scalar a(_grid);
    Scalar q1_cross_q2(_grid);
    Flux q2_cross_a(_grid);
    
    // TODO: Eventually, loop over several different values of q1, q2, and a
    q1 = 3.;
    q2 = 7.;
    a = 5.;
    q1_cross_q2 = CrossProduct( q1, q2 );
    q2_cross_a = CrossProduct( q2, a );
    double IPFlux = InnerProduct( q1, q2_cross_a );
    double IPScalar = InnerProduct( a, q1_cross_q2 );
    EXPECT_DOUBLE_EQ( IPFlux, IPScalar );
}

// ========
// = Curl =
// ========

// Test curl( a ) = 0 if a is constant in space
TEST_F(VectorOperationsTest, CurlOfConstantScalarIsZero) {
    Scalar a(_grid);
    a = 3;
    Flux q = Curl(a);
    EXPECT_ALL_X_EQ( 0, q(X,i,j) );
    EXPECT_ALL_Y_EQ( 0, q(Y,i,j) );
}

// Test curl( q ) = 0 if q is constant in space
TEST_F(VectorOperationsTest, CurlOfConstantFluxIsZero) {
    Flux q(_grid);
    q = 5;
    Scalar a = Curl(q);
    EXPECT_ALL_EQ( 0, a(i,j) );
}

// Test curl( (ax + c, by + d) ) = a - b (constant)
TEST_F(VectorOperationsTest, CurlOfLinearFlux) {
    double a = 7;
    double b = 11;
    double c = 13;
    double d = 19;
    Scalar u(_grid);
    Scalar v(_grid);
    Flux q(_grid);
    v = a * _x + c;
    u = b * _y + d;
    VelocityToFlux( u, v, q );
    Scalar CurlQ = Curl( q );
    EXPECT_ALL_INTERIOR_EQ( CurlQ(i,j), a - b );
}

// Test curl( ax + by + c ) = (b, -a)  (constant)
TEST_F(VectorOperationsTest, CurlOfLinearScalar) {
    double a = 7;
    double b = 11;
    double c = 13;
    Scalar f(_grid);
    f = a * _x + b * _y + c;
    Flux CurlF = Curl( f );

    double dx = _grid.getDx();
    EXPECT_ALL_X_EQ( CurlF(X,i,j), b * dx );
    EXPECT_ALL_Y_EQ( CurlF(Y,i,j), -a * dx );
}

// Test that < curl(q), a > = < q, curl(a) >
// for any Flux q, and any Scalar a with a = 0 on boundaries
TEST_F(VectorOperationsTest, CurlInnerProduct) {
    Flux q(_grid);
    Scalar a(_grid);

    // TODO: Eventually, loop over several different values of a and q,
    //       with a = 0 on boundaries
    InitializeSingleWavenumber( 1, 1, a );
    q = 3;
    Scalar curlQ = Curl( q );
    Flux curlA = Curl( a );
    double IPFlux = InnerProduct( q, curlA );
    double IPScalar = InnerProduct( curlQ, a );
    EXPECT_NEAR( IPFlux, IPScalar, 5e-15 );
}


/*--------------------------------------------------------------------------*/
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
	
	Scalar sf = Curl(fluxc);		
	EXPECT_ALL_EQ(sf(i,j), 0); 		
}

TEST_F(VectorOperationsTest, ScalarCurlOfFluxGivesZeroBCs) {		
	Scalar sf = Curl(_q);		
	for (int i=0; i<_nx+1; ++i) {
		EXPECT_DOUBLE_EQ( sf(i,0), 0 );
		EXPECT_DOUBLE_EQ( sf(i,_ny), 0 );
	}
	for (int j=0; j<_ny+1; ++j) {       
		EXPECT_DOUBLE_EQ( sf(0,j), 0 );
		EXPECT_DOUBLE_EQ( sf(_nx,j), 0 );  
	}                         
}	

// =================
// = Sin transform =
// =================

// Sine transform only counts inner grids and gives zero boundary conditions.
TEST_F(VectorOperationsTest, SineTransformGivesZeroBCs) {
	Scalar f_fft = SinTransform(_f);
	for (int i=0; i<_nx+1; ++i) {
		EXPECT_DOUBLE_EQ( f_fft(i,0), 0 );
		EXPECT_DOUBLE_EQ( f_fft(i,_ny), 0 );
	}
	for (int j=0; j<_ny+1; ++j) {       
		EXPECT_DOUBLE_EQ( f_fft(0,j), 0 );
		EXPECT_DOUBLE_EQ( f_fft(_nx,j), 0 );  
	} 
}

// Sine transform of zero gives zero.
TEST_F(VectorOperationsTest, SineTransformOfZeroEqualsZero) {	
	Scalar Szero(_grid);
	Szero = 0; 	
	
	Scalar S_fft = SinTransform(Szero);
	EXPECT_ALL_EQ( S_fft(i,j), 0 );
}

// Sine transform gives a Scalar of the same size as its argument
TEST_F(VectorOperationsTest, SineTransformOutputSize) {
    Scalar a(_grid);
    a = 0;
    
    Scalar b = SinTransform(a);
    
    // compare grids
    const Grid& g1 = a.getGrid();
    const Grid& g2 = b.getGrid();
    EXPECT_EQ( g1.getNy(), g2.getNy() );
    EXPECT_EQ( g1.getNy(), g2.getNy() );
}

// Do Sine transform twice and obtain the original data, upto the normalization parameter.
TEST_F(VectorOperationsTest, SineTransformTwiceGivesOriginalScalar) {
	Scalar f_fft = SinTransform(_f);
	Scalar f_f_fft = SinTransform(f_fft);
	double normalization = 2 * _nx * 2 * _ny; //normalization needed (see FFTW3 manual)
    double tolerance = 2e-14;
	for (int i = 1; i<_nx; ++i) {
		for (int j = 1; j<_ny; ++j) {
			EXPECT_NEAR( f_f_fft(i,j)/normalization, _f(i,j), tolerance ); 
		}	
	}		
}


// ================================
// = BoundaryVector inner product =
// ================================
TEST_F(VectorOperationsTest, BoundaryVectorInnerProduct) {
    const int n=10;
    BoundaryVector x(n);
    BoundaryVector y(n);
    BoundaryVector::index ind;
    
    x = 1;
    y = 1;
    EXPECT_DOUBLE_EQ( InnerProduct(x,y), 2.0*n );

    // loop over x-components
    for (ind = x.begin(X); ind != x.end(X); ++ind) {
        x(ind) = 2;
        y(ind) = 3;
    }
    // loop over y-components
    for (ind = x.begin(Y); ind != x.end(Y); ++ind) {
        x(ind) = 5;
        y(ind) = 7;
    }
    EXPECT_DOUBLE_EQ( InnerProduct(x,y), 6*n + 35*n );

}

#undef EXPECT_ALL_EQ
#undef EXPECT_ALL_X_EQ
#undef EXPECT_ALL_Y_EQ

} // namespace