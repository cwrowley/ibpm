#include <gtest/gtest.h>
#include "SingleWavenumber.h"
#include "VectorOperations.h"
#include "EllipticSolver.h"
#include "Array.h"
#include <math.h>

using namespace ibpm;
using Array::Array2;

namespace {
    
    const double tolerance = 1e-10;
    
#define EXPECT_ALL_EQ(a,b)                          \
    for (int lev=0; lev<ngrid; ++lev) {             \
        for (int i=1; i<nx; ++i) {                  \
            for (int j=1; j<ny; ++j) {              \
                EXPECT_NEAR( (a), (b), tolerance ); \
            }                                       \
        }                                           \
    }

//------------------------------------------------------------------------------
// Single-domain tests
//------------------------------------------------------------------------------
    
class EllipticSolverTest : public testing::Test {
protected:
    void TestPoisson( const Grid& grid );
    void TestHelmholtz( const Grid& grid, double alpha );
};

    
void EllipticSolverTest::TestPoisson( const Grid& grid ) {
    Scalar u( grid );
    Scalar f( grid );
    PoissonSolver poisson( grid );
    int nx = grid.Nx();
    int ny = grid.Ny();
    int ngrid = grid.Ngrid();

    for (int kx = 0; kx < nx; ++kx) {
        for (int ky = 0; ky < ny; ++ky) {
            InitializeSingleWavenumber( kx, ky, f );
            // Solve L u = f
            poisson.solve( f, u );
            Scalar Lu( grid );
            Laplacian( u, Lu );
            f.coarsify();
#ifdef DEBUG
            cout << "(kx, ky) = (" << kx << ", " << ky << ")" << endl;
            cout << "u:" << endl;
            u.print();
            cout << "Lu:" << endl;
            Lu.print();
            cout << "f:" << endl;
            f.print();
#endif
            EXPECT_ALL_EQ( f(lev,i,j), Lu(lev,i,j) );
        }
    }
}

void EllipticSolverTest::TestHelmholtz( const Grid& grid, double alpha ) {
    Scalar u( grid );
    Scalar f( grid );
    HelmholtzSolver helmholtz( grid, alpha );
    int nx = grid.Nx();
    int ny = grid.Ny();
    
    for (int kx = 0; kx < nx; ++kx) {
        for (int ky = 0; ky < ny; ++ky) {
            InitializeSingleWavenumber( kx, ky, f );
            // Solve (1 + alpha * Laplacian) u = f
            helmholtz.solve( f, u );
            Scalar Lu( grid );
            Laplacian( u, Lu );
            Lu *= alpha;
            Lu += u;
//            EXPECT_ALL_EQ( f(lev,i,j), Lu(lev,i,j) );
        }
    }
}

    TEST_F( EllipticSolverTest, PoissonSingleDomain ) {
        int nx = 4;
        int ny = 8;
        int ngrid = 1;
        double length = 0.4;
        Grid grid( nx, ny, ngrid, length, -1., -1. );
        TestPoisson( grid );
    }
    
    TEST_F( EllipticSolverTest, HelmholtzSingleDomain ) {
        int nx = 4;
        int ny = 8;
        int ngrid = 1;
        double length = 0.4;
        Grid grid( nx, ny, ngrid, length, -1., -1. );
        double alpha = 0.2;
        TestHelmholtz( grid, alpha );
    }
    
    TEST_F( EllipticSolverTest, PoissonMultiDomain ) {
        int nx = 4;
        int ny = 8;
        int ngrid = 2;
        double length = 0.4;
        Grid grid( nx, ny, ngrid, length, -1., -1. );
        TestPoisson( grid );
    }
    
    TEST_F( EllipticSolverTest, HelmholtzMultiDomain ) {
        int nx = 4;
        int ny = 8;
        int ngrid = 2;
        double length = 0.4;
        Grid grid( nx, ny, ngrid, length, -1., -1. );
        double alpha = 0.2;
        TestHelmholtz( grid, alpha );
    }
    
} // namespace
