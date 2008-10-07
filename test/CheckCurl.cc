// Test Curl operation

#include <iostream>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

void AdjointOfScalarCurl( const Flux& q, Scalar& f ) {
    // Loop over all Scalar fields
    for (int lev=0; lev<q.Ngrid(); ++lev) {
        for (int i=1; i<q.Nx(); ++i) {
            for (int j=1; j<q.Ny(); ++j) {
                Scalar e(q.getGrid()); // basis function
                e = 0;
                e(lev,i,j) = 1.;
                Flux curlE(q.getGrid()) ;
                Curl( e, curlE );
                f(lev,i,j) = InnerProduct( q, curlE ) / InnerProduct( e, e );
            }
        }
    }
}


int main() {
    bool printCurlScalar = true;
    bool printCurlFlux = true;
    int nx=8;
    int ny=8;
    int ngrid=2;
    double length=8; // dx = 1
    Grid grid( nx, ny, ngrid, length, -length/2, -length*ny/nx/2 );
    Scalar u(grid);
    Scalar v(grid);
    Flux q(grid);

    // Initialize scalars
    for (int lev=0; lev<ngrid; ++lev) {
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny; ++j) {
                u(lev,i,j) = grid.getXEdge(lev,i);
                v(lev,i,j) = grid.getYEdge(lev,j);
            }
        }
    }

    // Initialize flux:
    // velocity = ( -y, x )
    for (int lev=0; lev<ngrid; ++lev) {
        double dx = grid.Dx(lev);
        // X-flux
        for (int i=0; i<=nx; ++i) {
            for (int j=0; j<ny; ++j) {
                q(lev,X,i,j) = -2*q.y(lev,X,j)*q.y(lev,X,j) * dx;
            }
        }
        // Y-flux
        for (int i=0; i<nx; ++i) {
            for (int j=0; j<=ny; ++j) {
                q(lev,Y,i,j) = q.x(lev,Y,i)*q.x(lev,Y,i) * dx;
            }
        }
    }
    
    
    Flux curlU(grid);
    Flux curlV(grid);
    Curl( u, curlU );
    Curl( v, curlV );
    if (printCurlScalar) {
        cout << "u:" << endl;
        u.print();
        cout << "Curl(u):" << endl;
        curlU.print();
        
        cout << "v:" << endl;
        v.print();
        cout << "Curl(v):" << endl;
        curlV.print();
    }
    if (printCurlFlux) {
        cout << "q:" << endl;
        q.print();
        Scalar omega(grid);
        Curl( q, omega );
        cout << "Curl(q):" << endl;
        omega.print();
        Scalar omega2(grid);
        AdjointOfScalarCurl( q, omega2 );
        cout << "AdjointOfScalarCurl( q ):" << endl;
        omega2.print();
        Scalar err = omega - omega2;
        cout << "Error:" << endl;
        err.print();
    }    
    return 0;
}