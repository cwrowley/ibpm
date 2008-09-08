// SingleWavenumber.h - routine for initializing single-wavenumber fields
#include "Grid.h"
#include "Scalar.h"

using namespace ibpm;

// Set the output Scalar f equal to sin(kx * x) * sin(ky * y)
// with the specified wavenumbers in x and y
inline void InitializeSingleWavenumber(
    int xWavenumber,
    int yWavenumber,
    Scalar& f
    ) {
    const double pi = 4 * atan(1.);
    const int nx = f.getNx();
    const int ny = f.getNy();
    const double xLength = f.getXEdge(nx) - f.getXEdge(0);
    const double yLength = f.getYEdge(ny) - f.getYEdge(0);
    const double kx = xWavenumber * pi / xLength;
    const double ky = yWavenumber * pi / yLength;
    const double deltaX = f.getDx();
    
    for (int i=0; i <= nx; ++i) {
        double x = i * deltaX;
        for (int j=0; j <= ny; ++j) {
            double y = j * deltaX;
            f(i,j) = sin(kx * x) * sin(ky * y);
        }
    }
}
