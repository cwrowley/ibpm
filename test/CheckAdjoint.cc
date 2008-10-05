// Verify that VelocityToXFlux is adjoint of FluxToXVelocity

#include <iostream>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

void printY(Flux q) {
    for (int lev=0; lev<q.Ngrid(); ++lev) {
        cout << "level " << lev << endl;
        for (int j=q.Ny(); j>=0; --j) {
            for (int i=0; i<q.Nx(); ++i) {
                cout << q(lev,Y,i,j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

// Compute the adjoint of FluxToXVelocity, brute force:
// q(lev,i,j) = < q, e(lev,i,j) > / <e,e> where e is an orthogonal basis function
// x(lev,i,j) = < x, f(lev,i,j) > / <f,f>
// Mx (lev,i,j) = < Mx, e(lev,i,j) > = < x, M* e(lev,i,j) > / <e,e>
void AdjFluxToX( const Scalar& x, Flux& q ) {
    Flux e( q.getGrid() ); // basis function
    // loop over all Flux basis functions
    for (int lev=0; lev<q.Ngrid(); ++lev) {
        for (int i=0; i<=q.Nx(); ++i) {
            for (int j=0; j<q.Ny(); ++j) {
                e = 0;
                e(lev,X,i,j) = 1.;
                Scalar y(x.getGrid());
                FluxToXVelocity( e, y );
                double a = InnerProduct( x, y );
                double normsq = InnerProduct( e, e );
                q(lev,X,i,j) = a / normsq;
            }
        }
    }
}

// Compute the adjoint of XVelocityToFlux, brute force:
void AdjYVelToFlux( const Flux& q, Scalar& x ) {
    Scalar e( x.getGrid() ); // basis function
    // loop over all Scalar basis functions
    x = 0.;
    for (int lev=0; lev<x.Ngrid(); ++lev) {
        for (int i=1; i<x.Nx(); ++i) {
            for (int j=1; j<x.Ny(); ++j) {
                e = 0;
                e(lev,i,j) = 1.;
                Flux p(q.getGrid());
                p = 0;
                YVelocityToFlux( e, p );
                double a = InnerProduct( q, p );
                double normsq = InnerProduct( e, e );
                x(lev,i,j) = a / normsq;
            }
        }
    }
}

// In XVelocityToFlux, flux q(lev,i,j) depends on nearby scalar values
// Set x to the coefficient of each scalar value that influences q(lev,i,j)
void computeScalarDependencies(int lev, int i, int j, Scalar& x) {
    // Loop over all possible Scalars
    x = 0.;
    for (int ll=0; ll < x.Ngrid(); ++ll) {
        for (int ii=1; ii<x.Nx(); ++ii) {
            for (int jj=1; jj<x.Ny(); ++jj) {
                Scalar e( x.getGrid() );
                Flux q( x.getGrid() );
                e = 0.;
                e(ll,ii,jj) = 1.;
                XVelocityToFlux(e, q);
                x(ll,ii,jj) = q(lev,X,i,j);
            }
        }
    }
}


int main() {
    int nx=8;
    int ny=8;
    int ngrid=3;
    double length=0.1; // dx = 1
    Grid grid( nx, ny, ngrid, length, 0., 0. );
    Scalar err(grid);
    Flux err2(grid);
    Flux mag(grid);
    err2 = 0.;
    mag = 0.;
    
    // Loop over all possible Scalars
    for (int lev=0; lev<ngrid; ++lev) {
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny;++j) {
                Scalar x(grid);
                x = 0;
                x(lev,i,j) = 1;
                Flux q(grid);
                Flux q0(grid);
                q = 0;
                q0 = 0;
                AdjFluxToX( x, q0 );
                XVelocityToFlux( x, q );
                err(lev,i,j) = InnerProduct(q-q0,q-q0) / InnerProduct(q0,q0);
//                cout << "Scalar ("  << i << ", " << j << "):" << endl;
//                x.print();
//                cout << "Exact adjoint: " << endl;
//                printX(q0);
//                cout << "XVelocityToFlux:" << endl;
//                printX(q);
//                cout << "[return to continue]";
//                char c;
//                cin.get(c);
            }
        }
    }
    // Loop over all fluxes
    for (int lev=0; lev<ngrid; ++lev) {
        for (int ind=err2.begin(Y); ind < err2.end(Y); ++ind) {
            Flux e(grid);
            e = 0;
            e(lev,ind) = 1;
            Scalar x(grid);
            Scalar x0(grid);
            x = 0;
            x0 = 0;
            FluxToYVelocity( e, x0 );
            AdjYVelToFlux( e, x );
            mag(lev,ind) = InnerProduct(x0,x0);
            err2(lev,ind) = InnerProduct(x-x0,x-x0);
        }
    }
    cout << "Norm of FluxToYVelocity:" << endl;
    printY(mag);
    cout << "Errors in adjoint:" << endl;
    printY(err2);
    while (0) {
    //while (cin.good()) {
        cout << "Enter indices to look at influence coefs (lev,i,j): ";
        double lev, i, j;
        cin >> lev >> i >> j;
        Scalar x(grid);
        computeScalarDependencies(lev,i,j,x);
        cout << "Magnitude of FluxToVelocity:" << endl;
        printY(mag);
        cout << "Errors in adjoint:" << endl;
        printY(err2);
        err.print();
        cout << endl << "Influence coefs:" << endl;
        x.print();
        cout << "Error in adjoint: " << err(lev,i,j) << endl;
    }
    return 0;
}