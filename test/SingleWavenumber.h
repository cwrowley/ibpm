// SingleWavenumber.h - routine for initializing single-wavenumber fields
#include "Grid.h"
#include "Scalar.h"
#include "Array.h"
#include <math.h>
#include <iostream>

using namespace std;
using namespace ibpm;
using Array::Array2;

// Set the output Scalar f equal to sin(kx * x) * sin(ky * y)
// with the specified wavenumbers in x and y
inline void InitializeSingleWavenumber(
    int xWavenumber,
    int yWavenumber,
    Scalar& f
    ) {
    const double pi = 4 * atan(1.);
    const int nx = f.Nx();
    const int ny = f.Ny();
    const int ngrid = f.Ngrid();
    const double x0 = f.getXEdge(ngrid-1,0);
    const double y0 = f.getYEdge(ngrid-1,0);
    const double xLength = f.getXEdge(ngrid-1,nx) - x0;
    const double yLength = f.getYEdge(ngrid-1,ny) - y0;
    const double kx = xWavenumber * pi / xLength;
    const double ky = yWavenumber * pi / yLength;
    
    for (int lev=0; lev < ngrid; ++lev ) {
        for (int i=1; i < nx; ++i) {
            double x = f.getXEdge(lev,i);
            for (int j=1; j < ny; ++j) {
                double y = f.getYEdge(lev,j);
                f(lev,i,j) = sin(kx * (x-x0)) * sin(ky * (y-y0));
            }
        }
    }
}

// Set the output Array f equal to sin(kx * x) * sin(ky * y),
// where the domain is length pi in each direction
// The input array is assumed to have indices ( 1..nx, 1..ny )
inline void InitializeSingleWavenumber(
    int kx,
    int ky,
    Array2<double>& f
    ) {
    int nx = f.Nx() + 1;
    int ny = f.Ny() + 1;
    const double pi = 4 * atan(1.);
    double dx = pi / nx;
    double dy = pi / ny;
    for (int i=1; i < nx; ++i) {
        for (int j=1; j < ny; ++j) {
            f(i,j) = sin( kx * i * dx ) * sin( ky * j * dy );
        }
    }
}
