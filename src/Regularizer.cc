// Regularizer.cc
//
// Description:
// Implementation of the Regularizer class
//
// Author(s):
// Clancy Rowley
//
// Date: 7 Jul 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL: Header$

#include "Regularizer.h"
#include "Grid.h"
#include "Geometry.h"
#include "Flux.h"
#include <vector>
#include <math.h>

namespace ibpm {

Regularizer::Regularizer(const Grid& grid, const Geometry& geometry) :
    _grid(grid),
    _geometry(geometry)
{}

// number of cells over which delta function has support
const double deltaSupportRadius = 1.5;

// Return the value of the regularized delta function phi(r)
// where delta(x) \approx phi(x/h) / h, where h is the grid spacing
//
// From Roma, Peskin, and Berger, JCP 1999, eq (22)
// Note that the input r is given in numbers of cells
// (e.g. normalized by the cell width)
inline double deltaFunction(double r) {
    if (r > deltaSupportRadius) {
        return 0;
    }
    else {
        if (r <= 0.5)
            return (1 + sqrt(1 - 3*r*r)) / 3;
        else
            return (5 - 3*r - sqrt(1 - 3*(1-r)*(1-r)) ) / 6;
    }
}

// Update list of relationships between boundary points and cells, and the
// corresponding weights
// Checks only the finest grid level, level=0
void Regularizer::update() {
    Direction dir;
    Flux f(_grid);
    int i;
    Flux::index j;
    double dx, dy;
    double h = _grid.Dx();  // mesh spacing
    Association a;

    // Clear the list of associated Flux and BoundaryVector points
    _neighbors.clear();

    // Get the coordinates of the body
    BoundaryVector bodyCoords = _geometry.getPoints();

    // For each direction (x and y)
    for (dir = X; dir <= Y; ++dir) {
        // For each point on the boundary
        for (i = 0; i < bodyCoords.getNumPoints(); ++i) {
            // For each cell
            for (j = f.begin(dir); j != f.end(dir); ++j) {
                // Find x and y distances between boundary point and cell
                dx = fabs(f.x(0,j) - bodyCoords(X,i)) / h;
                dy = fabs(f.y(0,j) - bodyCoords(Y,i)) / h;
                // If cell is within the radius of support of delta function
                if ((dx < deltaSupportRadius) && (dy < deltaSupportRadius)) {
                    // Compute the weight factor
                    a.weight = deltaFunction(dx) * deltaFunction(dy);
                    a.fluxIndex = j;
                    a.boundaryIndex = bodyCoords.getIndex(dir,i);
                    // Add to list of associated cells
                    _neighbors.push_back(a);
                }
            }
        }
    }
}

Flux Regularizer::toFlux(const BoundaryVector& u1) const {
    // Allocate a new Flux field, initialized to zero
    Flux u2(_grid);
    u2 = 0;

    // For each association between cells and boundary points
    vector<Association>::const_iterator a;
    for (a = _neighbors.begin(); a != _neighbors.end(); ++a) {
        // add the weight factor times the boundary value to the flux
        u2(0,a->fluxIndex) += a->weight * u1(a->boundaryIndex);
    }

    // Multiply by grid spacing for correct dimension (vector -> Flux)
    u2 *= _grid.Dx();

    // Return the new flux field
    return u2;
}

BoundaryVector Regularizer::toBoundary(const Flux& u2) const {
    // Allocate a new BoundaryVector, initialized to zero
    BoundaryVector u1(_geometry.getNumPoints());
    u1 = 0;

    // For each association between cells and boundary points
    vector<Association>::const_iterator a;
    for (a = _neighbors.begin(); a != _neighbors.end(); ++a) {
        // Add the weight factor times the flux value to the boundary
        u1(a->boundaryIndex) += a->weight * u2(0,a->fluxIndex);
    }

    // Divide by grid spacing for correct dimension (Flux -> vector)
    u1 /= _grid.Dx();    

    // Return the new BoundaryVector
    return u1;
}

} // namespace

