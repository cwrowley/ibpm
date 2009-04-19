// CheckLaplacian -- driver routine to check Laplacian and Poisson solver

#include <iostream>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

PoissonSolver* solver;

bool equal( const Scalar& u, const Scalar& v ) {
    static double tol=1e-12;
    double err = 0;
    for (int lev=0; lev<u.Ngrid(); ++lev) {
        for (int i=1; i<u.Nx(); ++i) {
            for (int j=1; j<u.Ny(); ++j) {
                err += fabs( u(lev,i,j) - v(lev,i,j) );
            }
        }
    }
    return err < tol;
}

// Return true if PoissonSolver computes exact inverse of L(u),
// where u is zero at all points, except 1 at (lev,i,j)
void checkLinvLEqualsIdentity( const Grid& grid, int lev, int i, int j ) {
    Scalar u(grid);
    u = 0;
    u(lev,i,j) = 1;
    
    Scalar Lu = Laplacian( u );
    
    // Solve L w = Lu
    Scalar w(u.getGrid());
    solver->solve( Lu, w );
    
    // Check that w = u everywhere
    if ( ! equal(u, w) ) {
        cout << "Inverse not exact for (" << lev << ", " << i << ", " << j << ")" << endl;
//        cout << "u:" << endl;
//        u.print();
//        cout << "Lu:" << endl;
//        Lu.print();
//        cout << "Lu.coarsify()" << endl;
//        Lu.coarsify();
//        Lu.print();
//        cout << "w = Linv( Lu ):" << endl;
//        w.print();
//        cout << "Lw:" << endl;
//        Laplacian( w, Lu );
//        Lu.print();
    }
}

// Check L( Linv( u ) ) = u.coarsify()
void checkLofLinvEqualsCoarsify( const Grid& grid, int lev, int i, int j ) {
    Scalar f(grid);
    f = 0;
    f(lev,i,j) = 1;
    
    // Solve Lu = f
    Scalar u(grid);
    solver->solve( f, u );
    
    // Compute Lu
    Scalar Lu(grid);
    Laplacian( u, Lu );
    
//    f.coarsify();
    // Check that Lu = f everywhere
    if ( ! equal( Lu, f ) ) {
        cout << "Inverse not exact for (" << lev << ", " << i << ", " << j << ")" << endl;
    }
}

double InnerProduct2( const Scalar& f, const Scalar& g ) {
    double ip = 0;
    for (int lev=0; lev<f.Ngrid(); ++lev) {
        double dx2 = f.Dx(lev) * f.Dx(lev);
        for (int i=1; i<f.Nx(); ++i) {
            for (int j=1; j<f.Ny(); ++j) {
                ip += f(lev,i,j) * g(lev,i,j) * dx2;
            }
        }
    }
    return ip;
}

// Compute the adjoint of Laplacian, brute force:
// u(lev,i,j) = < u, e(lev,i,j) > / <e,e> where e is an orthogonal basis function
// Lu (lev,i,j) = < Lu, e(lev,i,j) > = < x, L* e(lev,i,j) > / <e,e>
void AdjLaplacian( const Scalar& f, Scalar& u ) {
    Scalar e( f.getGrid() ); // basis function
    // loop over all basis functions
    for (int lev=0; lev<f.Ngrid(); ++lev) {
        for (int i=1; i<f.Nx(); ++i) {
            for (int j=1; j<f.Ny(); ++j) {
                e = 0;
                e(lev,i,j) = 1.;
                Scalar y(f.getGrid());
                Laplacian( e, y );
                double a = InnerProduct2( f, y );
                double normsq = InnerProduct2( e, e );
                u(lev,i,j) = a / normsq;
            }
        }
    }
}

bool checkAdjLaplacian( const Grid& grid, int lev, int i, int j ) {
    Scalar u(grid);
    u = 0;
    u(lev,i,j) = 1;

    Scalar Lu(grid);
    Laplacian( u, Lu );
    Scalar LstarU(grid);
    AdjLaplacian( u, LstarU );
    
    if ( equal( Lu, LstarU ) ) {
        return true;
    }
    else {
//        cout << "Does not equal adjoint for (" << lev << ", " << i << ", " << j << ")" << endl;
//        cout << "u:" << endl;
//        u.print();
//        cout << "Lu:" << endl;
//        Lu.print();
//        cout << "L*u:" << endl;
//        LstarU.print();
        return false;
    }
}

// Return true if LCu = CLCu for u=0 everywhere except (lev,i,j)
bool checkLCEqualsCLC( const Grid& grid, int lev, int i, int j ) {
    Scalar u(grid);
    u = 0;
    u(lev,i,j) = 1;
    
    u.coarsify();
    Scalar LCu(grid);
    Laplacian( u, LCu );
    
    Scalar CLCu = LCu;
    CLCu.coarsify();
    return equal( CLCu, LCu );
}

int main() {
    int nx=8;
    int ny=8;
    int ngrid=2;
    double length=8; // dx = 1
    Grid grid( nx, ny, ngrid, length, -length/2, -length*ny/nx/2 );
    solver = new PoissonSolver( grid );
    
    Scalar status(grid);
    status = 2;
    
    for (int lev = 0; lev<ngrid; ++lev) {
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny; ++j) {
                // checkLinvLEqualsIdentity( grid, lev, i, j );
                // checkLofLinvEqualsCoarsify( grid, lev, i, j );
//                if ( checkAdjLaplacian( grid, lev, i, j ) ) {
                if ( checkLCEqualsCLC( grid, lev, i, j ) ) {
                    status(lev,i,j) = 0;
                }
                else {
                    status(lev,i,j) = 1;
                }
                cout << "." << flush;
            }
        }
    }
    cout << endl;
    status.print();

    // while (0) {
    while (cin.good()) {
        cout << "Enter indices to look at (lev,i,j): ";
        double lev, i, j;
        cin >> lev >> i >> j;

        Scalar u(grid);
        u = 0;
        u(lev,i,j) = 1;
        u.coarsify();
        Scalar LCu = Laplacian( u );
        Scalar CLCu = LCu;
        CLCu.coarsify();
                
        cout << "LCu:" << endl;
        LCu.print();
        
        cout << "CLCu:" << endl;
        CLCu.print();
        
        cout << endl;
    }
    
    
    return 0;
}