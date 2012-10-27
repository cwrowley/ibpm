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
#include <unistd.h>
#include <gtest/gtest.h>

using namespace ibpm;

namespace {

const double tolerance = 1e-10;

#define EXPECT_ALL_EQ(a,b)                      \
    for (int i=1; i<_nx; ++i) {                 \
        for (int j=1; j<_ny; ++j) {             \
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
        _ngrid(1),
        _timestep(0.123),
        _grid(_nx, _ny, _ngrid, 2., -1, -1) {

        double Reynolds = 100.;

        // Define two geometries: one with and one without a body
        RigidBody body;
        _emptyGeometry.addBody(body);

        int nPoints = 4;
        body.addLine_n( -0.75, 0, 0.75, 0, nPoints );
        _nonemptyGeometry.addBody(body);

        BaseFlow q0( _grid, 1.0, 0. );

        _modelWithNoBodies = new NavierStokesModel(
            _grid, _emptyGeometry, Reynolds, q0 );
        _modelWithNoBodies->init();
        
        _modelWithBodies = new NavierStokesModel(
            _grid, _nonemptyGeometry, Reynolds, q0 );
        _modelWithBodies->init();
    }

    ~ProjectionSolverTest() {
        delete _modelWithNoBodies;
        delete _modelWithBodies;
    }

    // Compute the linear term on the LHS of the projection equations:
    //   (1 - h/2 L) omega
    // where L is given by the associated NavierStokesModel
    Scalar ComputeLinearTerm(NavierStokesModel& model, const Scalar& omega) {
        Scalar result = Laplacian( omega );
        result *= model.getAlpha();
        result *= -_timestep / 2.;
        result += omega;
        return result;
    }
    
    // Verify that the projection equations are satisfied by the solution of
    // the given solver, for the given model:
    //  (1 - h/2 L) omega + h B f = a
    //                      C f = b
    void verify(
        NavierStokesModel& model,
        ProjectionSolver& solver
        ) {
        const int nPoints = model.getNumPoints();
        // Variables on rhs of projection equations
        Scalar a(_grid);
        InitializeSingleWavenumber( 1, 1, a );
        BoundaryVector b(nPoints);
        b = 3.;
        // Variables to solve for
        Scalar omega(_grid);
        omega = 0;
        BoundaryVector f(nPoints);
        f = 0;

        solver.solve( a, b, omega, f );

        // Verify equations are satisfied to specified tolerance
        Scalar lhs = ComputeLinearTerm( model, omega );
        Scalar forcingTerm( _grid );
        model.B( f, forcingTerm );
        lhs += _timestep * forcingTerm;
        EXPECT_ALL_EQ( a(0,i,j), lhs(0,i,j) );

        // Verify constraint is satisfied
        BoundaryVector constraint( nPoints );
        model.C( omega, constraint );
        EXPECT_ALL_BV_EQ( b(X,i), constraint(X,i), nPoints );
        EXPECT_ALL_BV_EQ( b(Y,i), constraint(Y,i), nPoints );
    }

    // data
    int _nx;
    int _ny;
    int _ngrid;
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
    ConjugateGradientSolver solver(_grid, *_modelWithNoBodies, _timestep, tolerance);
    verify( *_modelWithNoBodies, solver );
}

TEST_F(CGSolverTest, WithConstraints) {
    ConjugateGradientSolver solver( _grid, *_modelWithBodies, _timestep, tolerance);
    verify( *_modelWithBodies, solver );
}

TEST_F(CholeskySolverTest, NoConstraints) {
    CholeskySolver solver( _grid, *_modelWithNoBodies, _timestep );
    solver.init();
    verify( *_modelWithNoBodies, solver );
}

TEST_F(CholeskySolverTest, WithConstraints) {
    CholeskySolver solver( _grid, *_modelWithBodies, _timestep );
    solver.init();
    verify( *_modelWithBodies, solver );
}

TEST_F(CholeskySolverTest, SaveFile) {
    CholeskySolver solver( _grid, *_modelWithBodies, _timestep );
    solver.init();
    // Save a file
    bool success = solver.save("testSolver");
    EXPECT_EQ( true, success );

    // Load the file
    CholeskySolver newSolver( _grid, *_modelWithBodies, _timestep );
    success = newSolver.load("testSolver");
    EXPECT_EQ( true, success );

    verify( *_modelWithBodies, newSolver );
    
    // Attempt to load the file with a solver with the wrong number of points
    CholeskySolver differentBody( _grid, *_modelWithNoBodies, _timestep );
    success = differentBody.load("testSolver");
    EXPECT_EQ( false, success );
    
    // Attempt to load the file with a solver with the wrong timestep
    CholeskySolver differentTimestep( _grid, *_modelWithBodies, _timestep * 2. );
    success = differentTimestep.load("testSolver");
    EXPECT_EQ( false, success );
    
    unlink("testSolver.cholesky");
}

} // namespace
