#include <gtest/gtest.h>
#include "SingleWavenumber.h"
#include "VectorOperations.h"
#include "EllipticSolver2d.h"
#include "Array.h"
#include <math.h>
#include <iostream>

using namespace std;
using namespace ibpm;
using Array::Array2;

namespace {
    
const double tolerance = 1e-10;
    
#define EXPECT_ALL_EQ(a,b)                      \
    for (int i=1; i<_nx; ++i) {                 \
        for (int j=1; j<_ny; ++j) {             \
            EXPECT_NEAR( (a), (b), tolerance ); \
        }                                       \
    }
        
class EllipticSolver2dTest : public testing::Test {
protected:
    EllipticSolver2dTest() :
        _nx(4),
        _ny(8),
        _dx(0.1),
        _alpha(0.2),
        _poisson(_nx, _ny, _dx ),
        _helmholtz(_nx, _ny, _dx, _alpha)
    {}

    int _nx;
    int _ny;
    double _dx;
    double _alpha;
    PoissonSolver2d _poisson;
    HelmholtzSolver2d _helmholtz;
};
    
TEST_F( EllipticSolver2dTest, Poisson2d ) {
    Array2<double>  u( _nx-1, _ny-1, 1, 1 );
    Array2<double>  f( _nx-1, _ny-1, 1, 1 );
    Array2<double> Lu( _nx-1, _ny-1, 1, 1 );
    BC bc( _nx, _ny );
    bc = 0;
    
    for (int kx = 0; kx < _nx; ++kx) {
        for (int ky = 0; ky < _ny; ++ky) {
            InitializeSingleWavenumber( kx, ky, f );
            _poisson.solve( f, u );
            Laplacian( u, _dx, bc, Lu );
            EXPECT_ALL_EQ( f(i,j), Lu(i,j) );
        }
    }
}
    
TEST_F( EllipticSolver2dTest, Poisson2dConstWithBC ) {
    Array2<double> f(_nx-1, _ny-1, 1, 1);
    Array2<double> u(_nx-1, _ny-1, 1, 1);
    BC bc(_nx, _ny);
    // test f = 0, u = constant;
    f = 0.;
    const double val = 8.;
    bc = val;
    _poisson.solve( f, bc, u );
    EXPECT_ALL_EQ( val, u(i,j) );
}

// Test L u = 0 with bcs so that u = x + 2*y
TEST_F( EllipticSolver2dTest, Poisson2dLinearWithBC ) {
    Array2<double> f(_nx-1, _ny-1, 1, 1);
    Array2<double> u(_nx-1, _ny-1, 1, 1);
    Array2<double> u_exact(_nx-1, _ny-1, 1, 1);
    BC bc(_nx, _ny);

    f = 0.;
    
    // exact u (interior points)
    for (int i=1; i<_nx; ++i) {
        double x = i * _dx;
        for (int j=1; j<_ny; ++j) {
            double y = j * _dx;
            u_exact(i,j) = x + 2*y;
        }
    }
    
    // set boundary condition
    for (int i=0; i<=_nx; ++i) {
        double x = i * _dx;
        double y0 = 0;
        double y1 = _ny * _dx;
        bc.bottom(i) = x + 2*y0;
        bc.top(i) = x + 2*y1;
    }
    for (int j=0; j <= _ny; ++j) {
        double x0 = 0;
        double x1 = _nx * _dx;
        double y = j * _dx;
        bc.left(j) = x0 + 2*y;
        bc.right(j) = x1 + 2*y;
    }

    _poisson.solve( f, bc, u );
    EXPECT_ALL_EQ( u_exact(i,j), u(i,j) );
}

TEST_F( EllipticSolver2dTest, Helmholtz2d ) {
    Array2<double>  u( _nx-1, _ny-1, 1, 1 );
    Array2<double>  f( _nx-1, _ny-1, 1, 1 );
    Array2<double> Lu( _nx-1, _ny-1, 1, 1 );
    BC bc( _nx, _ny );
    bc = 0;
    
    for (int kx = 0; kx < _nx; ++kx) {
        for (int ky = 0; ky < _ny; ++ky) {
            InitializeSingleWavenumber( kx, ky, f );
            _helmholtz.solve( f, u );
            Laplacian( u, _dx, bc, Lu );
            Lu *= _alpha;
            Lu += u;
            EXPECT_ALL_EQ( f(i,j), Lu(i,j) );
        }
    }
}

// Test (1 + alpha * L) u = 0 with a constant (nonzero) boundary condition on u
TEST_F( EllipticSolver2dTest, Helmholtz2dWithConstBC ) {
    Array2<double>  f(_nx-1, _ny-1, 1, 1);
    Array2<double>  u(_nx-1, _ny-1, 1, 1);
    Array2<double> Lu(_nx-1, _ny-1, 1, 1);
    BC bc(_nx, _ny);
    // test f = 0, u = constant;
    f = 0.;
    const double val = 8.;
    bc = val;
    _helmholtz.solve( f, bc, u );
    Laplacian( u, _dx, bc, Lu );
    EXPECT_ALL_EQ( u(i,j), -_alpha * Lu(i,j) );
}

// Test (1 + alpha * L) u = 0 with a non-constant boundary condition on u
TEST_F( EllipticSolver2dTest, Helmholtz2dWithLinearBC ) {
    Array2<double>  f(_nx-1, _ny-1, 1, 1);
    Array2<double>  u(_nx-1, _ny-1, 1, 1);
    Array2<double> Lu(_nx-1, _ny-1, 1, 1);
    BC bc(_nx, _ny);
    
    f = 0.;
    
    // set boundary condition
    for (int i=0; i<=_nx; ++i) {
        double x = i * _dx;
        double y0 = 0;
        double y1 = _ny * _dx;
        bc.bottom(i) = x + 2*y0;
        bc.top(i) = x + 2*y1;
    }
    for (int j=0; j <= _ny; ++j) {
        double x0 = 0;
        double x1 = _nx * _dx;
        double y = j * _dx;
        bc.left(j) = x0 + 2*y;
        bc.right(j) = x1 + 2*y;
    }
    
    _helmholtz.solve( f, bc, u );
    Laplacian( u, _dx, bc, Lu );
    EXPECT_ALL_EQ( u(i,j), -_alpha * Lu(i,j) );
}

} // namespace
