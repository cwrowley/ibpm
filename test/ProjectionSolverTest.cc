#include "Grid.h"
#include "Flux.h"
#include "Scalar.h"
#include "NavierStokesModel.h"
#include "RigidBody.h"
#include "Geometry.h"
#include "ProjectionSolver.h"
#include "ConjugateGradientSolver.h"
#include <gtest/gtest.h>

namespace {

class ProjectionSolverTest : public testing::Test {
protected:
    ProjectionSolverTest() : 
        _nx(10),
        _ny(20),
        _timestep(1.),
        _grid(_nx, _ny, 2, -1, -2),
        _q0(_grid) {

        // Choose Reynolds number such that linear term is Laplacian
        _Reynolds = 1. / ( _grid.getDx() * _grid.getDx() );

        RigidBody body;
        body.addPoint(0,0);
        _geom.addBody(body);
        _q0 = Flux::UniformFlow( _grid, 1.0, 0. );
        _model = new NonlinearNavierStokes( _grid, _geom, _Reynolds, _q0 );
        _solver = new ConjugateGradientSolver( *_model, _timestep, 1e-10 );
    }

    // data
    int _nx;
    int _ny;
    double _timestep;
    double _Reynolds;
    Grid _grid;
    Flux _q0;
    Geometry _geom;
    NavierStokesModel* _model;
    ProjectionSolver* _solver;
};

TEST_F(ProjectionSolverTest, NoConstraints) {

    // Define a solver for an empty Geometry: no boundary points
    Geometry emptyGeom;
    NonlinearNavierStokes model( _grid, emptyGeom, _Reynolds, _q0 );
    ConjugateGradientSolver solver( model, _timestep, 1e-10 );

    BoundaryVector emptyVector(0);
    Scalar rhs(_grid);
    rhs = 1;
    Scalar gamma(_grid);
    gamma = 0;
    BoundaryVector f(0);
    
    solver.solve( rhs, emptyVector, gamma, f );

    // TODO: Actually check that rhs is correct
    // Note: the code below relies on the protected method Ainv -- bad
    
    // Scalar gamma2 = solver.Ainv( rhs );
    // for (int i=0; i <= _nx; ++i) {
    //     for (int j=0; j <= _ny; ++j) {
    //         EXPECT_DOUBLE_EQ( gamma2(i,j), gamma(i,j) );
    //     }
    // }
}

TEST_F(ProjectionSolverTest, True) {
    EXPECT_DOUBLE_EQ(1., 1.);
}

} // namespace
