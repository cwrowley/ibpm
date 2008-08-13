#include "Grid.h"
#include "Scalar.h"
#include "Flux.h"
#include "BoundaryVector.h"
#include "VectorOperations.h"
#include <gtest/gtest.h>

namespace {

class VectorOperationsTest : public testing::Test {
protected:
    VectorOperationsTest() : 
        _nx(4),
        _ny(8),
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

/*----------------------------------------------------------------------------*/
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

// Sine transform only counts inner grids and gives zero boundary conditions.
TEST_F(VectorOperationsTest, SineTransformGivesZeroBCs) {
	Scalar f_fft = sinTransform(_f);
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
	
	Scalar S_fft = sinTransform(Szero);
	EXPECT_ALL_EQUAL(S_fft(i,j), 0);
}

// Do Sine transform twice and obtain the original data, upto the normalization parameter.
TEST_F(VectorOperationsTest, SineTransformTwiceGivesOriginalScalar) {
	Scalar f_fft = sinTransform(_f);
	Scalar f_f_fft = sinTransform(f_fft);
	double normalization = 2 * _nx * 2 * _ny; //normalization needed (see FFTW3 manual)
	for (int i = 1; i<_nx; ++i) {
		for (int j = 1; j<_ny; ++j) {
			EXPECT_DOUBLE_EQ( f_f_fft(i,j)/normalization, _f(i,j) ); 
		}	
	}		
}

//TEST_F(VectorOperationsTest, SineTransformUsingDirectDST) {
//	Scalar f_fft = sinTransform(_f);
//	// to do: Compare the fft (RODFT00) result with the result by using DST-I definition...  
//}

//TEST_F(VectorOperationsTest, SineTransformOfSineFunction) {
// should give a const..  
//}

// X direction average of a cosnt Flux object gives a constant Scalar, except at boundary.
TEST_F(VectorOperationsTest, ScalarFluxXAverageOfConstFlux) {
	Flux fluxc(_grid);
	double cx = 8;
	double cy = 3;
	fluxc = cx;
	
	for (int i=0; i<_nx; ++i) {
		for (int j=0; j<_ny+1; ++j) {
			fluxc(Y,i,j) = cy;
		}
	}
	
	Scalar sc = fluxXAverage(fluxc); 
	for (int i=0; i<_nx+1; ++i) {
		for (int j=1; j<_ny; ++j){
			EXPECT_DOUBLE_EQ( sc(i,j), cx); 
		}
		EXPECT_DOUBLE_EQ( sc(i,0), 0.5 * cx);
		EXPECT_DOUBLE_EQ( sc(i,_ny), 0.5 * cx);
	}	
}

// X direction average of a simple linear Flux object.
TEST_F(VectorOperationsTest, ScalarFluxXAverageOfSimpleFlux) {
	Scalar sc = fluxXAverage(_p); // p.x = 2i+3j
	for (int i=0; i<_nx+1; ++i) {
		for (int j=1; j<_ny; ++j){
			EXPECT_DOUBLE_EQ( sc(i,j), 2 * i + 3 * j - 1.5); 
		}
		EXPECT_DOUBLE_EQ( sc(i,0), i);
		EXPECT_DOUBLE_EQ( sc(i,_ny), 0.5 * (2 * i + 3 * (_ny-1)));
	}	
}

// Y direction average of a cosnt Flux object gives a Scalar, 
// constant everywhere except at boundary.
TEST_F(VectorOperationsTest, ScalarFluxYAverageOfConstFlux) {
	Flux fluxc(_grid);
	double cx = 8;
	double cy = 3;
	fluxc = cx;
	
	for (int i=0; i<_nx; ++i) {
		for (int j=0; j<_ny+1; ++j) {
			fluxc(Y,i,j) = cy;
		}
	}
	
	Scalar sc = fluxYAverage(fluxc); 
	for (int j = 0; j < _ny+1; ++j) {
		for (int i= 1; i<_nx; ++i) {
			EXPECT_DOUBLE_EQ( sc(i,j), cy); 
		}
		EXPECT_DOUBLE_EQ( sc(0,j), 0.5 * cy);
		EXPECT_DOUBLE_EQ( sc(_nx,0), 0.5 * cy);
	}	
}

// Y direction average of a simple linear Flux object.
TEST_F(VectorOperationsTest, ScalarFluxYAverageOfSimpleFlux) {
	Scalar sc = fluxYAverage(_p); // p.y = i*j-10
	for (int j = 0; j < _ny+1; ++j) {
		for (int i= 1; i<_nx; ++i) {
			EXPECT_DOUBLE_EQ( sc(i,j), i * j - 10 - 0.5 * j); 
		}
		EXPECT_DOUBLE_EQ( sc(0, j), -5);
		EXPECT_DOUBLE_EQ( sc(_nx,j), 0.5 * (_nx - 1) * j -5);
	}	
}

//  Scalar generated by the cross product of two costant fluxes _q and _p is 
//  constant everywhere (except at boundary).
TEST_F(VectorOperationsTest, ScalarCrossProductOfTwoConstantFluxes) {
	Flux fluxc1(_grid);
	double c1x = 8;
	double c1y = 3;
	fluxc1 = c1x;
	for (int i=0; i<_nx; ++i) {
		for (int j=0; j<_ny+1; ++j) {
			fluxc1(Y,i,j) = c1y;
		}
	}
	
	Flux fluxc2(_grid);
	double c2x = 7;
	double c2y = 9;
	fluxc2 = c2x;
	for (int i=0; i<_nx; ++i) {
		for (int j=0; j<_ny+1; ++j) {
			fluxc2(Y,i,j) = c2y;
		}
	}
	
	Scalar f = crossproduct(fluxc1, fluxc2);
	double a =  c1x * c2y - c1y * c2x; // cross product
	for (int i = 1; i<_nx; ++i){
		for (int j = 1; j<_ny; ++j){
			EXPECT_DOUBLE_EQ( f(i,j), a); // at inner nodes
		}
		EXPECT_DOUBLE_EQ( f(i,0), 0.5 * a);   // at bottom 		
		EXPECT_DOUBLE_EQ( f(i,_ny), 0.5 * a);  // at top
	}
	
	for (int j = 1; j<_ny; ++j) {		
		EXPECT_DOUBLE_EQ( f(0,j), 0.5 * a);  // at left		
		EXPECT_DOUBLE_EQ( f(_nx,j), 0.5 * a);	// at right 	
	}
	
	// at four corners:
	EXPECT_DOUBLE_EQ( f(0, 0), 0.5 * 0.5 * a);
	EXPECT_DOUBLE_EQ( f(0, _ny), 0.5 * 0.5 * a);
	EXPECT_DOUBLE_EQ( f(_nx, 0), 0.5 * 0.5 * a);
	EXPECT_DOUBLE_EQ( f(_nx, _ny), 0.5 * 0.5 * a);
}


//  Scalar generated by the cross product of two simple fluxes _q and _p
// (with the zero 'ghost fluxes' assumption).
TEST_F(VectorOperationsTest, ScalarCrossProductOfTwoSimpleFluxes) {
	Scalar f = crossproduct(_q, _p);
	// q.x = 3i*nx+ 5(j+1)i; q.y=4ij+ny+cos(i)
	// p.x = 2i+3j; py=ij-10.
	double a;
	for (int i = 1; i<_nx; ++i){
		for (int j = 1; j<_ny; ++j){
			// cross product of two fluxes' averages:
			a = (3 * i * _nx + 5 * i * j + 2.5 * i) * (i * j - 10 - 0.5 * j)
				- (4 * i * j + _ny + 0.5 * (cos(i) + cos(i-1)) -2 * j) * (2 * i + 3 * j - 1.5);
			EXPECT_DOUBLE_EQ( f(i,j), a); // at inner nodes
		}
		int j = 0; // at bottom  
		a = 0.5 * (3 * i * _nx + 5 * i ) * (i * j - 10 - 0.5 * j)
			- (4 * i * j + _ny + 0.5 * (cos(i) + cos(i-1)) -2 * j) * i;
		EXPECT_DOUBLE_EQ( f(i,j), a);   
		j = _ny; // at top  
		a = 0.5 * (3 * i * _nx + 5 * _ny * i ) * (i * j - 10 - 0.5 * j)
			- (4 * i * j + _ny + 0.5 * (cos(i) + cos(i-1)) -2 * j) * 0.5 * (2 * i + 3 * (_ny - 1));
		EXPECT_DOUBLE_EQ( f(i,j), a);  
	}
	
	for (int j = 1; j<_ny; ++j) {
		int i = 0; // at left 
		a = (3 * i * _nx + 5 * i * j + 2.5 * i) * (-5)
			- 0.5 * (_ny + cos(0)) * (2 * i + 3 * j - 1.5);
		EXPECT_DOUBLE_EQ( f(i,j),  a);  
		i = _nx; // at right
		a = (3 * i * _nx + 5 * i * j + 2.5 * i) *  (0.5 * (_nx - 1) * j - 5)
			- 0.5 * (4 * (_nx - 1) * j + _ny + cos(_nx - 1)) * (2 * i + 3 * j - 1.5);
		EXPECT_DOUBLE_EQ( f(i,j),  a);	 	
	}
	
	// at four corners:
	int i = 0; int j = 0; 
	EXPECT_DOUBLE_EQ( f(i, j),  
		0.5 * (3 * i * _nx + 5 * i ) * (-5) - 0.5 * (_ny + cos(0)) * i);
	i = 0; j = _ny;
	EXPECT_DOUBLE_EQ( f(i, j),  
		0.5 * (3 * i * _nx + 5 * _ny * i ) * (-5) 
			- 0.5 * (_ny + cos(0)) * 0.5 * (2 * i + 3 * (_ny - 1)));
	i = _nx; j = 0;
	EXPECT_DOUBLE_EQ( f(i, j),  
		0.5 * (3 * i * _nx + 5 * i ) * (0.5 * (_nx - 1) * j - 5) 
			- 0.5 * (4 * (_nx - 1) * j + _ny + cos(_nx - 1)) * i);
	i = _nx; j = _ny;
	EXPECT_DOUBLE_EQ( f(i, j),  
		0.5 * (3 * i * _nx + 5 * _ny * i ) * (0.5 * (_nx - 1) * j - 5) 
			- 0.5 * (4 * (_nx - 1) * j + _ny + cos(_nx - 1)) * 0.5 * (2 * i + 3 * (_ny - 1)));	
}

/*----------------------------------------------------------------------------*/
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

// Cross product of a Flux and a Scalar gives zero BCs for the resulting flux.
TEST_F(VectorOperationsTest, FluxCrossProductZeroBCs) {
	Flux fluxf = crossproduct(_q, _g);	
	for (int j=0; j<_ny; ++j) {
		EXPECT_DOUBLE_EQ( fluxf(X,0,j), 0 );
		EXPECT_DOUBLE_EQ( fluxf(X,_nx,j), 0 );
	}
	for (int i=0; i<_nx; ++i) {       
		EXPECT_DOUBLE_EQ( fluxf(Y,i,0), 0 );
		EXPECT_DOUBLE_EQ( fluxf(Y,i,_ny), 0 );  
	} 	
}

// Cross product of a constant Flux and a constant Scalar
TEST_F(VectorOperationsTest, FluxCrossProductOfConstFluxAndScalar) {
	Scalar scalarc(_grid);
	double cs = 10;
	scalarc = cs;
	Flux fluxc(_grid);
	double cx = 8;
	double cy = 3;
	fluxc = cx;
	for (int i=0; i<_nx; ++i) {
		for (int j=0; j<_ny+1; ++j) {
			fluxc(Y,i,j) = cy;
		}
	}
	Flux fluxf = crossproduct(fluxc, scalarc);
	double a;
	for (int i=1; i<_nx; ++i) {
		for (int j=0; j<_ny; ++j){
			a = cy * cs;
			EXPECT_DOUBLE_EQ( fluxf(X,i,j), a); 
		}
	}
	for (int i=0; i<_nx; ++i) {
		for (int j=1; j<_ny; ++j){
			a = -cx * cs;
			EXPECT_DOUBLE_EQ( fluxf(Y,i,j), a); 
		}
	}
}

// Cross product of a simple Flux and a simple Scalar.
TEST_F(VectorOperationsTest, FluxCrossProductOfSimpleFluxAndScalar) {
	Flux fluxf = crossproduct(_p, _g); // _p.X = 2i+3j; _p.Y=ij-10; _g=-5i^2+3ij+2j^2.
	double a;
	for (int i=1; i<_nx; ++i) {
		for (int j=0; j<_ny; ++j){
			a = 0.5*((i*j-10-j*0.5)*g(i,j)+(i*(j+1)-10-(j+1)*0.5)*g(i,j+1));
			EXPECT_DOUBLE_EQ( fluxf(X,i,j), a); 
		}
	}
	for (int i=0; i<_nx; ++i) {
		for (int j=1; j<_ny; ++j){
			a = -0.5*((2*i+3*j-3*0.5)*g(i,j)+(2*(i+1)+3*j-3*0.5)*g(i+1,j));
			EXPECT_DOUBLE_EQ( fluxf(Y,i,j), a); 
		}
	}
}

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

#undef EXPECT_ALL_EQUAL
#undef EXPECT_ALL_X_EQUAL
#undef EXPECT_ALL_Y_EQUAL

} // namespace