#include "Array.h"
#include "BC.h"
#include "Grid.h"
#include "Scalar.h"
#include "Flux.h"
#include "BoundaryVector.h"
#include "VectorOperations.h"
#include "SingleWavenumber.h"
#include <gtest/gtest.h>

using Array::Array2;

using namespace ibpm;

namespace {

class VectorOperationsTest : public testing::Test {
protected:
    VectorOperationsTest() : 
        _nx(8),
        _ny(8),
        _ngrid(3),
        _grid(_nx, _ny, _ngrid, 2, -1, -3),
        _f(_grid),
        _g(_grid),
        _x(_grid),
        _y(_grid),
        _p(_grid),
        _q(_grid) {
            
        // Initialize Scalars _f, _g, _x, _y
        for (int lev = 0; lev < _ngrid; ++lev) {
            for (int i=1; i<_nx; ++i) {
                for (int j=1; j<_ny; ++j) {
                    _f(lev,i,j) = f(lev,i,j);
                    _g(lev,i,j) = g(lev,i,j);       
                    _x(lev,i,j) = _grid.getXEdge(lev,i);
                    _y(lev,i,j) = _grid.getYEdge(lev,j);
                }
            }
        }
        
        // Initialize Flux _q
        for (int lev = 0; lev < _ngrid; ++lev) {
            for (int i=0; i<_nx+1; ++i) {
                for (int j=0; j<_ny; ++j) {
                    _q(lev,X,i,j) = qx(lev,i,j);
                    _p(lev,X,i,j) = px(lev,i,j);                
                }
            }
        }
        for (int lev = 0; lev < _ngrid; ++lev) {
            for (int i=0; i<_nx; ++i) {
                for (int j=0; j<_ny+1; ++j) {
                    _q(lev,Y,i,j) = qy(lev,i,j);
                    _p(lev,Y,i,j) = py(lev,i,j);
                }
            }
        }
    }
    
    // functions to test
    // TODO: Pick better functions, with analytic expressions (in x and y)
    // (Single wavenumber for the scalar fields?)
    inline double f(int lev, int i, int j) {
        return 0.5 * i * i * _nx + 2 * j * (i + 1) + cos((double)j) * (_ny + 1);
    }

    inline double g(int lev, int i, int j) {
        return -5 * i * i + 3 * i * j + 2 * j * j;
    }
        
    inline double qx(int lev, int i, int j) {
        return 3 * i * _nx + 5 * (j +1 )* i;
    }
    
    inline double qy(int lev, int i, int j) {
        return 4 * i *  j   + _ny + cos((double)i);
    }
    
    inline double px(int lev, int i, int j) {
        return 2*i + 3*j;
    }

    inline double py(int lev, int i, int j) {
        return i*j - 10;
    }
    
    // Area of domain included by inner product of X-fluxes
    // (right and left edges excluded)
    double fluxXArea() {
        double Lx = _grid.getXCenter(_ngrid-1,_nx-1) - _grid.getXCenter(_ngrid-1,0);
        double Ly = _grid.getYEdge(_ngrid-1,_ny) - _grid.getYEdge(_ngrid-1,0);
        return Lx * Ly / ( _grid.Dx() * _grid.Dx() );
    }

    // Area of domain included by inner product of Y-flues
    // (top and bottom edges excluded)
    double fluxYArea() {
        double Lx = _grid.getXEdge(_ngrid-1,_nx) - _grid.getXEdge(_ngrid-1,0);
        double Ly = _grid.getYCenter(_ngrid-1,_ny-1) - _grid.getYCenter(_ngrid-1,0);
        return Lx * Ly / ( _grid.Dx() * _grid.Dx() );
    }
    
    double scalarArea() {
        double Lx = _grid.getXCenter(_ngrid-1,_nx-1) - _grid.getXCenter(_ngrid-1,0);
        double Ly = _grid.getYCenter(_ngrid-1,_ny-1) - _grid.getYCenter(_ngrid-1,0);
        return Lx * Ly;
    }
    
    // data
    int _nx;
    int _ny;
    int _ngrid;
    Grid _grid;
    Scalar _f;
    Scalar _g;
    Scalar _x;
    Scalar _y;
    Flux _p;
    Flux _q;
};

// for Scalars
#define EXPECT_ALL_EQ(a,b)                      \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for (int i=1; i<_nx; ++i) {             \
            for (int j=1; j<_ny; ++j) {         \
                EXPECT_DOUBLE_EQ( (a), (b) );   \
            }                                   \
        }                                       \
    }

// for Fluxes
#define EXPECT_ALL_X_EQ(a,b)                    \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for ( int i=0; i<_nx+1; ++i ) {         \
            for ( int j=0; j<_ny; ++j ) {       \
                EXPECT_DOUBLE_EQ( (a), (b) );   \
            }                                   \
        }                                       \
    }

// for Fluxes
#define EXPECT_ALL_Y_EQ(a,b)                    \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for ( int i=0; i<_nx; ++i ) {           \
            for ( int j=0; j<_ny+1; ++j ) {     \
                EXPECT_DOUBLE_EQ( (a), (b) );   \
            }                                   \
        }                                       \
    }

TEST_F(VectorOperationsTest, ScalarDotProductSymmetric) {
    EXPECT_DOUBLE_EQ(InnerProduct(_f,_g), InnerProduct(_f,_g) );        
}

TEST_F(VectorOperationsTest, ScalarDotProductDomainArea) {
    Scalar one( _grid );
    one = 1.;
    EXPECT_DOUBLE_EQ( scalarArea(), InnerProduct( one, one ) );
}

TEST_F(VectorOperationsTest, FluxDotProductSymmetric) {
    EXPECT_DOUBLE_EQ( InnerProduct(_q,_p), InnerProduct(_p, _q) );
}

TEST_F(VectorOperationsTest, ScalarDotProductZeroVector) {
    Scalar one( _grid );
    one = 1.;
    Scalar zero( _grid );
    zero = 0.;
    double ip = InnerProduct( one, zero );
    EXPECT_DOUBLE_EQ( 0., ip );
    ip = InnerProduct( zero, one );
    EXPECT_DOUBLE_EQ( 0., ip );
}

TEST_F(VectorOperationsTest, FluxDotProductZeroVector) {
    Flux one( _grid );
    one = 1.;
    Flux zero( _grid );
    zero = 0.;
    double ip = InnerProduct( one, zero );
    EXPECT_DOUBLE_EQ( 0., ip );
    ip = InnerProduct( zero, one );
    EXPECT_DOUBLE_EQ( 0., ip );
}

TEST_F(VectorOperationsTest, FluxDotProductDomainArea) {
    Flux h( _grid );
    Flux l( _grid );
    for (int lev=0; lev<_ngrid; ++lev) {
        double dx = exp2(lev);
        for (int i=0; i<_nx+1; ++i) {
            for (int j=0; j<_ny; ++j) {
                h(lev,X,i,j) = 1 * dx;
                l(lev,X,i,j) = 1 * dx;
            }
        }
        for (int i=0; i<_nx; ++i) {
            for (int j=0; j<_ny+1; ++j) {
                h(lev,Y,i,j) = 0 * dx;
                l(lev,Y,i,j) = 0 * dx;
            }
        }
    }
    EXPECT_DOUBLE_EQ( InnerProduct(h,l), fluxXArea() );

    for (int lev=0; lev<_ngrid; ++lev) {
        double dx = exp2(lev);
        for (int i=0; i<_nx+1; ++i) {
            for (int j=0; j<_ny; ++j) {
                h(lev,X,i,j) = 0 * dx;
                l(lev,X,i,j) = 0 * dx;
            }
        }
        for (int i=0; i<_nx; ++i) {
            for (int j=0; j<_ny+1; ++j) {
                h(lev,Y,i,j) = 1 * dx;
                l(lev,Y,i,j) = 1 * dx;
            }
        }
    }
    EXPECT_DOUBLE_EQ( InnerProduct(h,l), fluxYArea() );
}

// =============================
// = Flux to velocity and back =
// =============================

// Set X and Y Flux values at the coarsest grid level to xval and yval, resp.
void SetFluxBoundary(double xval, double yval, Flux& q) {
    int lev = q.Ngrid()-1;
    // X flux
    for (int i=0; i<=q.Nx(); ++i) {
        q(lev,X,i,0) = xval;
        q(lev,X,i,q.Ny()-1) = xval;
    }
    for (int j=1; j<q.Ny()-1; ++j) {
        q(lev,X,0,j) = xval;
        q(lev,X,q.Nx(),j) = xval;
    }
    
    // Y flux
    for (int j=0; j<=q.Ny(); ++j) {
        q(lev,Y,0,j) = yval;
        q(lev,Y,q.Nx()-1,j) = yval;
    }
    for (int i=1; i<q.Nx()-1; ++i) {
        q(lev,Y,i,0) = yval;
        q(lev,Y,i,q.Ny()) = yval;
    }
}

// If u = const, then corresponding X-Flux should be constant too
TEST_F(VectorOperationsTest, ConstantVelocityToFlux) {
    Scalar u(_grid);
    Scalar v(_grid);
    Flux q(_grid);
    double xval = 4.;
    double yval = 12.;
    double dx = _grid.Dx();
    q = 0.;
    u = xval;
    v = yval;
    
    XVelocityToFlux( u, q );
    YVelocityToFlux( v, q );
    // fudge boundaries
    double dx_coarse = dx * exp2(_ngrid-1);
    SetFluxBoundary( xval * dx_coarse, yval * dx_coarse, q );
#ifdef DEBUG
    cout << "q:" << endl;
    q.print();
#endif
    EXPECT_ALL_X_EQ( xval * dx * exp2(lev), q(lev,X,i,j) );
    EXPECT_ALL_Y_EQ( yval * dx * exp2(lev), q(lev,Y,i,j) );
}

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
    
    q = 5.;
    // TODO: Eventually, test several different values of q
    Scalar a = CrossProduct( q, q );
    EXPECT_ALL_EQ( 0, a(lev,i,j) );
}

// Test that (u,0) x (0,v) = u v
TEST_F(VectorOperationsTest, CrossProductOfOrthogonalVectors) {
    Scalar u(_grid);
    Scalar v(_grid);
    Flux q1(_grid);
    Flux q2(_grid);

    // TODO: Eventually, test several different values of u, v
    u = 3;
    v = 5;
    q1 = 0;
    q2 = 0;
    XVelocityToFlux( u, q1 );
    YVelocityToFlux( v, q2 );
    Scalar cross = CrossProduct( q1, q2 );
    Scalar err = u * v - cross;
    // set boundary elements to zero
    for (int lev=0; lev<err.Ngrid(); ++lev) {
        for (int i=1; i<err.Nx(); ++i) {
            err(lev,i,1) = 0;
            err(lev,i,err.Ny()-1) = 0;
        }
        for (int j=1; j<u.Ny(); ++j) {
            err(lev,1,j) = 0;
            err(lev,u.Nx()-1,j) = 0;
        }
    }
    double normsq = InnerProduct( err, err );
#ifdef DEBUG
    cout << "u:" << endl;
    u.print();
    cout << "v:" << endl;
    v.print();
    cout << "u x v:" << endl;
    cross.print();
    cout << "err:" << endl;
    err.print();
#endif
    EXPECT_NEAR( 0., normsq, 1e-28 );
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
    EXPECT_ALL_EQ( q1q2(lev,i,j), -q2q1(lev,i,j) );
}

// Test that < a, q1 x q2 > = < q1, q2 x a >
// for all Scalars a, Fluxes q1, q2
TEST_F(VectorOperationsTest, CrossProductTripleProduct) {
    Flux q1(_grid);
    Flux q2(_grid);
    Scalar a(_grid);
    Scalar q1_cross_q2(_grid);
    Flux q2_cross_a(_grid);

    q1 = _q;
    q2 = _p;
    a = 5.;
    q1_cross_q2 = CrossProduct( q1, q2 );
    q2_cross_a = CrossProduct( q2, a );
#ifdef DEBUG
    cout << "q1:" << endl;
    q1.print();
    cout << "q2:" << endl;
    q2.print();
    cout << "a:" << endl;
    a.print();
    cout << "q2 x a:" << endl;
    q2_cross_a.print();
    cout << "q1 x q2:" << endl;
    q1_cross_q2.print();
#endif
    double IPFlux = InnerProduct( q1, q2_cross_a );
    double IPScalar = InnerProduct( a, q1_cross_q2 );
    double tol = 1e-14 * IPFlux;
    EXPECT_NEAR( IPFlux, IPScalar, tol );
}

// ========
// = Curl =
// ========

// Test curl( a ) = 0 if a is constant in space
TEST_F(VectorOperationsTest, CurlOfConstantScalarIsZero) {
    Scalar a(_grid);
    a = 3;
    Flux q = Curl(a);
    // fudge outer boundary
    SetFluxBoundary( 0, 0, q );
    EXPECT_ALL_X_EQ( 0, q(lev,X,i,j) );
    EXPECT_ALL_Y_EQ( 0, q(lev,Y,i,j) );    
}

// Test curl( q ) = 0 if q is constant in space
TEST_F(VectorOperationsTest, CurlOfConstantFluxIsZero) {
    Flux q(_grid);
    q = 5;
    Scalar a = Curl(q);
    EXPECT_ALL_EQ( 0, a(lev,i,j) );
}

// Set boundary values of f to specified value val
void SetScalarBoundary( double val, Scalar& f ) {
    int lev = f.Ngrid() - 1;
    for (int i=1; i<f.Nx(); ++i) {
        f(lev,i,1) = val;
        f(lev,i,f.Ny()-1) = val;
    }
    for (int j=1; j<f.Ny(); ++j) {
        f(lev,1,j) = val;
        f(lev,f.Nx()-1,j) = val;
    }
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
    // fudge outer boundaries
    SetScalarBoundary( a - b, CurlQ );
#ifdef DEBUG
    q.print();
    CurlQ.print();
#endif
    EXPECT_ALL_EQ( a - b, CurlQ(lev,i,j) );
}

// Test curl( ax + by + c ) = (b, -a)  (constant)
TEST_F(VectorOperationsTest, CurlOfLinearScalar) {
    double a = 7;
    double b = 11;
    double c = 13;
    Scalar f(_grid);
    f = a * _x + b * _y + c;
    Flux CurlF = Curl( f );
    
    // fudge outer boundary of coarsest grid:
    // should not be correct value, since Curl assumes zero boundary conditions
    // on coarsest grid
    // So set to correct value here
    double dx = _grid.Dx() * exp2(_ngrid-1);
    SetFluxBoundary( b * dx, -a * dx, CurlF );
    
    dx = _grid.Dx();
#ifdef DEBUG
    f.print();
    CurlF.print();
    cout << "X: " << b * dx << endl;
    cout << "Y: " << -a * dx << endl;
#endif
    EXPECT_ALL_X_EQ( CurlF(lev,X,i,j), b * dx * exp2(lev) );
    EXPECT_ALL_Y_EQ( CurlF(lev,Y,i,j), -a * dx * exp2(lev) );
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

TEST_F(VectorOperationsTest, ScalarCurlOfConstantEqualsZero) {      
    // Test 1: Curl of a non-zero constant Flux object
    // (w/ different Flux.x and Flux.y).
    Flux fluxc(_grid);
    double cx = 8;
    double cy = 5;
    fluxc = cx;
    for (int lev=0; lev<_ngrid; ++lev) {
        for (int i=0; i<_nx; ++i) {
            for (int j=0; j<_ny+1; ++j) {
                fluxc(lev,Y,i,j) = cy;
            }
        }
    }
    
    Scalar sf = Curl(fluxc);
    EXPECT_ALL_EQ( 0., sf(lev,i,j) );      
}

TEST_F(VectorOperationsTest, LaplacianOfLinearEqualsZero) {
    // construct linear field u = x + 2*y
    Array2<double> u(_nx-1, _ny-1, 1, 1);
    Array2<double> Lu(_nx-1, _ny-1, 1, 1);
    BC bc(_nx,_ny);

    // interior points
    for (int i=1; i<_nx; ++i) {
        double x = _grid.getXEdge(0,i);
        for (int j=1; j<_ny; ++j) {
            double y = _grid.getYEdge(0,j);
            u(i,j) = x + 2*y;
        }
    }

    // boundary points
    for (int i=0; i<=_nx; ++i) {
        double x = _grid.getXEdge(0,i);
        double y0 = _grid.getYEdge(0,0);
        double y1 = _grid.getYEdge(0,_ny);
        bc.bottom(i) = x + 2*y0;
        bc.top(i) = x + 2*y1;
    }
    for (int j=0; j <= _ny; ++j) {
        double x0 = _grid.getXEdge(0,0);
        double x1 = _grid.getXEdge(0,_nx);
        double y = _grid.getYEdge(0,j);
        bc.left(j) = x0 + 2*y;
        bc.right(j) = x1 + 2*y;
    }
    Laplacian( u, _grid.Dx(), bc, Lu );
    EXPECT_ALL_EQ( 0., Lu(i,j) );
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
