#include "Grid.h"
#include "Flux.h"
#include "Scalar.h"
#include "NavierStokesModel.h"
#include "RigidBody.h"
#include "Geometry.h"
#include "ProjectionSolver.h"
#include "ConjugateGradientSolver.h"
#include "SingleWavenumber.h"
#include <gtest/gtest.h>

namespace {

const double tolerance = 1e-10;
    
class ProjectionSolverTest : public testing::Test {
protected:
    ProjectionSolverTest() : 
        _nx(6),
        _ny(6),
        _timestep(0.123),
        _grid(_nx, _ny, 2., -1, -1),
        _q0(_grid) {

        // Choose Reynolds number such that linear term is Laplacian
        _Reynolds = 1.;

        RigidBody body;
        body.addPoint(0.25,0.25);
        _geom.addBody(body);
        _nPoints = _geom.getNumPoints();
        // _q0 = Flux::UniformFlow( _grid, 1.0, 0. );
        _q0 = 0.;
        _model = new NonlinearNavierStokes( _grid, _geom, _Reynolds, _q0 );
        _solver = new ConjugateGradientSolver( *_model, _timestep, tolerance);
    }

    // Compute the linear term on the LHS of the projection equations:
    //   (1 + h/2 L) gamma
    // where L is given by the associated NavierStokesModel
    Scalar ComputeLinearTerm(NavierStokesModel* model, const Scalar& gamma) {
        Scalar result = model->S( gamma );
        result *= *( model->getLambda() );
        result = model->Sinv( result );
        result *= -_timestep / 2.;
        result += gamma;
        return result;
    }
    
    // data
    int _nx;
    int _ny;
    int _nPoints;
    double _timestep;
    double _Reynolds;
    Grid _grid;
    Flux _q0;
    Geometry _geom;
    NavierStokesModel* _model;
    ProjectionSolver* _solver;
};

#define EXPECT_ALL_EQ(a,b)                      \
	for (int i=0; i<_nx+1; ++i) {               \
		for (int j=0; j<_ny+1; ++j) {           \
			EXPECT_NEAR( (a), (b), tolerance ); \
		}                                       \
	}

#define EXPECT_ALL_BV_EQ(a,b)               \
    for (int i=0; i<_nPoints; ++i) {        \
        EXPECT_NEAR( (a), (b), tolerance ); \
    }

// TEST_F(ProjectionSolverTest, AOfAinvEqualsIdentity ) {
//     Scalar gamma(_grid);
//     InitializeSingleWavenumber( 1, 1, gamma );
//     
//     Scalar ATimesGamma = ComputeLinearTerm( _model, gamma );
//     Scalar AinvAGamma = _solver->Ainv( ATimesGamma );
//     EXPECT_ALL_EQ( AinvAGamma(i,j), gamma(i,j) );
// }

TEST_F(ProjectionSolverTest, NoConstraints) {

    // Define a solver for an empty Geometry: no boundary points
    Geometry emptyGeom;
    NonlinearNavierStokes model( _grid, emptyGeom, _Reynolds, _q0 );
    ConjugateGradientSolver solver( model, _timestep, 1e-10 );

    BoundaryVector emptyVector(0);
    Scalar rhs(_grid);
    InitializeSingleWavenumber( 1, 1, rhs );
    Scalar gamma(_grid);
    gamma = 0;
    BoundaryVector f(0);
    
    solver.solve( rhs, emptyVector, gamma, f );
    Scalar lhs = ComputeLinearTerm( &model, gamma );
    EXPECT_ALL_EQ( rhs(i,j), lhs(i,j) );

}

TEST_F(ProjectionSolverTest, WithConstraints) {
    // Variables on rhs of projection equations
    Scalar a(_grid);
    InitializeSingleWavenumber( 1, 1, a );
    BoundaryVector b(_geom.getNumPoints());
    b = 3.;
    // Variables to solve for
    Scalar gamma(_grid);
    gamma = 0;
    BoundaryVector f(_geom.getNumPoints());
    f = 0;

    _solver->solve( a, b, gamma, f );

    // Verify equations are satisfied to specified tolerance
    Scalar lhs = ComputeLinearTerm( _model, gamma );
    Scalar forcingTerm = _model->B( f );
    lhs += _timestep * forcingTerm;
    EXPECT_ALL_EQ( a(i,j), lhs(i,j) );

    // Verify constraint is satisfied
    BoundaryVector constraint = _model->C( gamma );
    EXPECT_ALL_BV_EQ( b(X,i), constraint(X,i) );
    EXPECT_ALL_BV_EQ( b(Y,i), constraint(Y,i) );
    
}

} // namespace
