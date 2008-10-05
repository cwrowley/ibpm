#include <iostream>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

void printX(Flux q) {
    for (int lev=0; lev<q.Ngrid(); ++lev) {
        cout << "level " << lev << endl;
        for (int j=q.Ny()-1; j>=0; --j) {
            for (int i=0; i<=q.Nx(); ++i) {
                cout << q(lev,X,i,j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

// In FluxToXVelocity, scalar x(lev,i,j) depends on nearby flux values
// Set q to the coefficient of each flux value that influences x(lev,i,j)
void computeDependencies(int lev, int i, int j, Flux& q) {
    // Loop over all possible X-fluxes
    q = 0.;
    for (int l=0; l < q.Ngrid(); ++l) {
        for (int ind=q.begin(X); ind != q.end(X); ++ind) {
            Flux e( q.getGrid() );
            Scalar x( q.getGrid() );
            e = 0.;
            e(l,ind) = 1.;
            FluxToXVelocity(e, x);
            q(l,ind) = x(lev,i,j);
        }
    }
}

int main() {
    int nx=8;
    int ny=8;
    int ngrid=3;
    double length=8; // dx = 0.1
    Grid grid( nx, ny, ngrid, length, 0., 0. );
    
    // Loop over all possible X-fluxes
    if (0) {
    for (int lev=0; lev<ngrid; ++lev) {
        for (int i=0; i<=nx; ++i) {
            for (int j=0; j<ny;++j) {
                Flux q(grid);
                Scalar x(grid);
                q = 0;
                q(lev,X,i,j) = 1;
                cout << "Flux (" << i << ", " << j << "):" << endl;
                printX(q);
                cout << "Scalar:" << endl;
                FluxToXVelocity( q, x );
                x.print();
                cout << "[return to continue]";
                char c;
                cin.get(c);
            }
        }
    }
    }
    
    // See what a constant field looks like
    Scalar u(grid);
    u = 1.;
    cout << "u:" << endl;
    u.print();
    Flux p(grid);
    XVelocityToFlux( u, p );
    cout << "ToFlux(u):" << endl;
    FluxToXVelocity( p, u );
    printX( p );
    cout << "ToScalar(ToFlux(u)):" << endl;
    u.print();
    
    for (int lev=0; lev<ngrid; ++lev) {
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny; ++j) {
                Flux q(grid);
                computeDependencies(lev,i,j,q);
                cout << "Scalar (" << lev << ", " << i << ", " << j << ") depends on Fluxes:"<<endl;
                printX(q);
                cout << "[return to continue]" << endl;
                char c;
                cin.get(c);
            }
        }
    }
                
    return 0;
}