#include "Grid.h"
#include "Flux.h"
#include "Scalar.h"
#include "NavierStokesModel.h"
#include "RigidBody.h"
#include "Geometry.h"
#include "ProjectionSolver.h"
#include "ConjugateGradientSolver.h"
#include "CholeskySolver.h"
#include "SingleWavenumber.h"
#include <gtest/gtest.h>

using namespace ibpm;

namespace {

const double tolerance = 1e-10;
    
Scalar ComputeLinearTerm(NavierStokesModel& model, const Scalar& gamma);

#define EXPECT_ALL_EQ(a,b)                      \
    for (int i=0; i<_nx+1; ++i) {               \
        for (int j=0; j<_ny+1; ++j) {           \
            EXPECT_NEAR( (a), (b), tolerance ); \
        }                                       \
    }

#define EXPECT_ALL_BV_EQ(a,b,nPoints)               \
    for (int i=0; i < nPoints; ++i) {        \
        EXPECT_NEAR( (a), (b), tolerance ); \
    }

class ProjectionSolverTest : public testing::Test {
protected:
    ProjectionSolverTest() :
        _nx(6),
        _ny(6),
        _timestep(0.123),
        _grid(_nx, _ny, 2., -1, -1) {

        // Choose Reynolds number such that linear term is Laplacian
        double Reynolds = 1.;

        // Define two geometries: one with and one without a body
        RigidBody body;
        _emptyGeometry.addBody(body);

        int nPoints = 4;
        body.addLine( -0.75, 0, 0.75, 0, nPoints );
        _nonemptyGeometry.addBody(body);

        Flux q0 = Flux::UniformFlow( _grid, 1.0, 0. );

        _modelWithNoBodies = new NonlinearNavierStokes(
            _grid, _emptyGeometry, Reynolds, q0 );
        _modelWithNoBodies->init();
        
        _modelWithBodies = new NonlinearNavierStokes(
            _grid, _nonemptyGeometry, Reynolds, q0 );
        _modelWithBodies->init();
    }

    ~ProjectionSolverTest() {
        delete _modelWithNoBodies;
        delete _modelWithBodies;
    }

    // Compute the linear term on the LHS of the projection equations:
    //   (1 - h/2 L) gamma
    // where L is given by the associated NavierStokesModel
    Scalar ComputeLinearTerm(NavierStokesModel& model, const Scalar& gamma) {
        Scalar result = model.S( gamma );
        result *= model.getLambda();
        result = model.Sinv( result );
        result *= -_timestep / 2.;
        result += gamma;
        return result;
    }
    
    // Verify that the projection equations are satisfied by the solution of
    // the given solver, for the given model:
    //  (1 - h/2 L) gamma + h B f = a
    //                      C f = b
    void verify(
        NavierStokesModel& model,
        ProjectionSolver& solver
        ) {
        const int nPoints = model.getGeometry().getNumPoints();
        // Variables on rhs of projection equations
        Scalar a(_grid);
        InitializeSingleWavenumber( 1, 1, a );
        BoundaryVector b(nPoints);
        b = 3.;
        // Variables to solve for
        Scalar gamma(_grid);
        gamma = 0;
        BoundaryVector f(nPoints);
        f = 0;

        solver.solve( a, b, gamma, f );

        // Verify equations are satisfied to specified tolerance
        Scalar lhs = ComputeLinearTerm( model, gamma );
        Scalar forcingTerm = model.B( f );
        lhs += _timestep * forcingTerm;
        EXPECT_ALL_EQ( a(i,j), lhs(i,j) );

        // Verify constraint is satisfied
        BoundaryVector constraint = model.C( gamma );
        EXPECT_ALL_BV_EQ( b(X,i), constraint(X,i), nPoints );
        EXPECT_ALL_BV_EQ( b(Y,i), constraint(Y,i), nPoints );
    }

    // data
    int _nx;
    int _ny;
    double _timestep;
    Grid _grid;
    Geometry _emptyGeometry;
    Geometry _nonemptyGeometry;
    NavierStokesModel* _modelWithNoBodies;
    NavierStokesModel* _modelWithBodies;
};

typedef ProjectionSolverTest CGSolverTest;
typedef ProjectionSolverTest CholeskySolverTest;

TEST_F(CGSolverTest, NoConstraints) {
    ConjugateGradientSolver solver(*_modelWithNoBodies, _timestep, tolerance);
    verify( *_modelWithNoBodies, solver );
}

TEST_F(CGSolverTest, WithConstraints) {
    ConjugateGradientSolver solver( *_modelWithBodies, _timestep, tolerance);
    verify( *_modelWithBodies, solver );
}

TEST_F(CholeskySolverTest, NoConstraints) {
    CholeskySolver solver( *_modelWithNoBodies, _timestep );
    solver.init();
    verify( *_modelWithNoBodies, solver );
}

TEST_F(CholeskySolverTest, WithConstraints) {
    CholeskySolver solver( *_modelWithBodies, _timestep );
    solver.init();
    verify( *_modelWithBodies, solver );
}

TEST_F(CholeskySolverTest, SaveFile) {
    CholeskySolver solver( *_modelWithBodies, _timestep );
    solver.init();
    bool success = solver.save("test.cholesky");
    EXPECT_EQ( true, success );

    CholeskySolver newSolver( *_modelWithBodies, _timestep );
    success = newSolver.load("test.cholesky");
    EXPECT_EQ( true, success );

    verify( *_modelWithBodies, newSolver );
}



} // namespace
