#include "RigidBody.h"
#include "Geometry.h"
#include "NavierStokesModel.h"
#include "SingleWavenumber.h"
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
        BaseFlow q0( _grid, magnitude, angle );
        _model = new NavierStokesModel( _grid, _geom, _Reynolds, q0 );
        _model->init();
    }

    Scalar linearTerm(const Scalar& g) {
        Scalar a( g.getGrid() );
        Laplacian( g, a );
        a *= _model->getAlpha();
        return a;
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

#define EXPECT_ALL_EQ(a,b)                  \
    for (int i=1; i<_nx; ++i) {             \
        for (int j=1; j<_ny; ++j) {         \
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
/*

// If the base flow is zero and
//    q = computeFlux(omega)
// then curl(q) should equal omega
TEST_F( NavierStokesModelTest, GammaToFlux ) {
    Scalar omega(_grid);
    Flux q(_grid);
    Scalar curlQ(_grid);

    // create a model with no base flow
    NavierStokesModel model( _grid, _geom, _Reynolds );
    model.init();
    // For a range of wavenumbers
    for (int xWavenumber = 0; xWavenumber < _nx; ++xWavenumber ) {
        for (int yWavenumber = 0; yWavenumber < _ny; ++yWavenumber) {
            InitializeSingleWavenumber( xWavenumber, yWavenumber, omega );
            model.computeFlux( omega, q );
            Curl( q, curlQ );
            EXPECT_ALL_EQ( omega(0,i,j), curlQ(0,i,j) );
            
            // Now put into a state vector and check refreshState()
            NavierStokesModel* modelp = &model;
            State x( _grid, _geom.getNumPoints() );
            x.q = 0.;
            x.omega = omega;
            modelp->refreshState( x );
            EXPECT_ALL_X_EQ( q(0,X,i,j), x.q(0,X,i,j) );
            EXPECT_ALL_Y_EQ( q(0,Y,i,j), x.q(0,Y,i,j) );
        }
    }
}
*/
} // namespace
