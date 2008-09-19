#include "RigidBody.h"
#include "Geometry.h"
#include "NavierStokesModel.h"
#include <gtest/gtest.h>
#include <iostream>

using namespace std;
using namespace ibpm;

namespace {

// Print a Scalar to standard out
void print(const Scalar& f);

class NavierStokesModelTest : public testing::Test {
protected:
    NavierStokesModelTest() :
        _nx(4),
        _ny(4),
        _ngrid(1),
        _length(4.*atan(1.)), // length is pi
        _xOffset(-1),
        _yOffset(-3),
        _grid( _nx, _ny, _ngrid, _length, _xOffset, _yOffset ) {
        
        // Choose Reynolds number such that linear term is Laplacian
        _Reynolds = 1.;
        // Add the geometry
        RigidBody body;
        body.addPoint(0,0);
        _geom.addBody(body);
        double magnitude = 1.;
        double angle = 0.;
        Flux q0 = Flux::UniformFlow( _grid, magnitude, angle );
        _model = new NonlinearNavierStokes( _grid, _geom, _Reynolds, q0 );
        _model->init();
    }

    Scalar linearTerm(const Scalar& g) {
        Scalar a( g.getGrid() );
        Laplacian( g, a );
        a *= _model->getAlpha();
        return a;
    }
    
    // Return Laplacian of the given Scalar, returning zero at the boundary
    Scalar testLaplacian(const Scalar& g) {
        Scalar Lg(g.getGrid());
        Lg = 0;
        for (int i=1; i<_nx; ++i) {
            for (int j=1; j<_ny; ++j) {
                Lg(i,j) = g(i-1,j) + g(i+1,j) + g(i,j-1) + g(i,j+1) -4*g(i,j);
                Lg(i,j) /= _grid.Dx() * _grid.Dx();
            }
        }
        return Lg;
    }

    // data
    int _nx;
    int _ny;
    int _ngrid;
    double _length;
    double _xOffset;
    double _yOffset;
    const Grid _grid;
    Geometry _geom;
    double _Reynolds;
    NavierStokesModel* _model;
};

// Set the output Scalar f equal to sin(kx * x) * sin(ky * y)
// with the specified wavenumbers in x and y
void InitializeSingleWavenumber(
    int xWavenumber,
    int yWavenumber,
    Scalar& f
    ) {
    const double pi = 4 * atan(1.);
    const Grid& grid = f.getGrid();
    const int nx = grid.Nx();
    const int ny = grid.Ny();
    const double xLength = grid.getXEdge(nx) - grid.getXEdge(0);
    const double yLength = grid.getYEdge(ny) - grid.getYEdge(0);
    const double kx = xWavenumber * pi / xLength;
    const double ky = yWavenumber * pi / yLength;
    const double deltaX = grid.Dx();
    
    for (int i=0; i <= nx; ++i) {
        double x = i * deltaX;
        for (int j=0; j <= ny; ++j) {
            double y = j * deltaX;
            f(i,j) = sin(kx * x) * sin(ky * y);
        }
    }
}

#define EXPECT_ALL_EQ(a,b)                  \
    for (int i=0; i<_nx+1; ++i) {           \
        for (int j=0; j<_ny+1; ++j) {       \
            EXPECT_NEAR((a), (b), 3e-15);   \
        }                                   \
    }

#define EXPECT_ALL_X_EQ(a,b)              \
    for ( int i=0; i<_nx+1; ++i ) {       \
        for ( int j=0; j<_ny; ++j ) {     \
            EXPECT_DOUBLE_EQ( (a), (b) ); \
        }                                 \
    }
    
#define EXPECT_ALL_Y_EQ(a,b)              \
    for ( int i=0; i<_nx; ++i ) {         \
        for ( int j=0; j<_ny+1; ++j ) {   \
            EXPECT_DOUBLE_EQ( (a), (b) ); \
        }                                 \
    }
    
//TEST_F( NavierStokesModelTest, SOfSinvEqualsIdentity ) {
//    Scalar gamma(_grid);
//    InitializeSingleWavenumber( 1, 1, gamma );
//    Scalar Sgamma = _model->S( gamma );
//    Scalar SinvSgamma = _model->Sinv( Sgamma );
//    EXPECT_ALL_EQ( SinvSgamma(i,j), gamma(i,j) );
//#ifdef DEBUG
//    cout << "gamma:" << endl;
//    gamma.print();
//    cout << "Sgamma:" << endl;
//    Sgamma.print();
//    cout << "SinvSgamma:" << endl;
//    SinvSgamma.print();
//#endif
//
//}

// If Reynolds number = 1, linear term should be Laplacian
TEST_F( NavierStokesModelTest, LinearTerm ) {
    // Test Laplacian of sinusoids
    Scalar gamma(_grid);
    
    // For a range of wavenumbers
    for (int xWavenumber = 0; xWavenumber < _nx; ++xWavenumber ) {
        for (int yWavenumber = 0; yWavenumber < _ny; ++yWavenumber) {
            InitializeSingleWavenumber( xWavenumber, yWavenumber, gamma );
            Scalar L_gamma = linearTerm( gamma );
            Scalar LaplacianGamma = testLaplacian( gamma );
            EXPECT_ALL_EQ( LaplacianGamma(i,j), L_gamma(i,j) );            
        }
    }
}

// If the base flow is zero and
//    q = computeFlux(gamma)
// then curl(q) should equal omega == gamma / dx^2
TEST_F( NavierStokesModelTest, GammaToFlux ) {
    Scalar gamma(_grid);
    Flux q(_grid);
    Scalar curlQ(_grid);
    double dx2 = _grid.Dx() * _grid.Dx();

    // create a model with no base flow
    NonlinearNavierStokes model( _grid, _geom, _Reynolds );
    model.init();
    // For a range of wavenumbers
    for (int xWavenumber = 0; xWavenumber < _nx; ++xWavenumber ) {
        for (int yWavenumber = 0; yWavenumber < _ny; ++yWavenumber) {
            InitializeSingleWavenumber( xWavenumber, yWavenumber, gamma );
            model.computeFlux( gamma, q );
            Curl( q, curlQ );
            EXPECT_ALL_EQ( gamma(i,j), curlQ(i,j) * dx2 );
            
            // Now put into a state vector and check refreshState()
            Model* modelp = &model;
            State x( _grid, _geom.getNumPoints() );
            x.q = 0.;
            x.gamma = gamma;
            modelp->refreshState( x );
            EXPECT_ALL_X_EQ( q(X,i,j), x.q(X,i,j) );
            EXPECT_ALL_Y_EQ( q(Y,i,j), x.q(Y,i,j) );
        }
    }
}

} // namespace
