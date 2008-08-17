#include "RigidBody.h"
#include "Geometry.h"
#include "NavierStokesModel.h"
#include <gtest/gtest.h>

namespace {

class NavierStokesModelTest : public testing::Test {
protected:
    NavierStokesModelTest() :
        _nx(10),
        _ny(20),
        _length(2),
        _xOffset(-1),
        _yOffset(-3),
        _grid( _nx, _ny, _length, _xOffset, _yOffset ) {
        
        // Choose Reynolds number such that linear term is Laplacian
        _Reynolds = 1. / ( _grid.getDx() * _grid.getDx() );

        // Add the geometry
        // NOTE: For now, Geometry stub returns a single point at the origin
        RigidBody body;
        body.addPoint(0,0);
        _geom.addBody(body);
        double magnitude = 1.;
        double angle = 0.;
        Flux q0 = Flux::UniformFlow( _grid, magnitude, angle );
        _model = new NonlinearNavierStokes( _grid, _geom, _Reynolds, q0 );        
    }

    Scalar linearTerm(const Scalar& g) {
        Scalar a = _model->S( g );
        a *= *( _model->getLambda() );
        a = _model->Sinv( a );
        return a;
    }

    // data
    int _nx;
    int _ny;
    double _length;
    double _xOffset;
    double _yOffset;
    const Grid _grid;
    Geometry _geom;
    double _Reynolds;
    NavierStokesModel* _model;
};

#define EXPECT_ALL_EQ(a,b)                  \
	for (int i=0; i<_nx+1; ++i) {           \
		for (int j=0; j<_ny+1; ++j) {       \
			EXPECT_DOUBLE_EQ((a), (b));     \
		}                                   \
	}

// If Reynolds number = 1, linear term should be Laplacian
TEST_F( NavierStokesModelTest, LinearTerm ) {
    // Test Laplacian of sinusoids
    Scalar gamma(_grid);
    Scalar L_gammaExact(_grid);
    const double pi = 4 * atan(1.);
    
    // For a range of wavenumbers
    for (int wavenumber = 0; wavenumber < 2; ++wavenumber ) {
        double length = _grid.getXEdge(_nx);
        double k = pi * wavenumber / length;
        
        // Define a Scalar as sin(kx)
        for (int i=0; i < _nx+1; ++i) {
            for (int j=0; j < _ny+1; ++j) {
                double x = _grid.getXEdge(i);
                gamma(i,j) = sin( k * x );
                L_gammaExact(i,j) = - k * k * gamma(i,j);
            }
        }
        
        // Take its Laplacian
        // Scalar L_gamma = linearTerm( gamma );
        Flux f = curl(gamma);
        Scalar L_gamma = curl(f);

        // Compare them
        EXPECT_ALL_EQ( L_gammaExact(i,j), L_gamma(i,j) );
    }
}

} // namespace
