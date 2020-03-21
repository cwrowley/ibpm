#ifndef _REGULARIZER_H_
#define _REGULARIZER_H_

#include <vector>
#include "Flux.h"
#include "BoundaryVector.h"

using std::vector;

namespace ibpm {

class Grid;
class Geometry;

/*!
\file Regularizer.h
\class Regularizer

\brief Define regularization operations between grid and boundary data

Regularization (smearing) lifts values defined on the boundary to values
defined on a grid, and interpolation takes data defined on a grid and
interpolates it to the boundary.

Specifically, if \f$ u(x) \f$ are values of u defined on a grid, and if
\f$ u(\xi) \f$ denotes values on a boundary, then regularization defines

\f[
u(x) = \int_\Omega u(\xi) \delta(x-\xi)\,d\xi
\f]

where a discrete approximation of the delta function is used, and
interpolation defines

\f[
u(\xi) = \int_\Omega u(x) \delta(x-\xi)\,dx
\f]

These integrals are discretized using a regularized version of the delta
function, with finite support, as in (14) of Taira & Colonius (J Comput Phys,
2007).

\author Clancy Rowley
\author $LastChangedBy$
\date  7 Jul 2008
\date $LastChangedDate$
\version $Revision$
*/

class Regularizer {
public:
    /// Constructor
    Regularizer(const Grid& grid, const Geometry& geometry);
    
    /// Destructor
    ~Regularizer() {}
    
    /// Update operators, for instance when the position of the bodies changes
    void update();
    
    /// \brief Smear boundary data to grid.
    /// In particular, if u1 denotes the vectors along the boundary,
    /// and u2 denotes the velocity vectors in the 2d domain, computes
    /// a discrete approximation to
    ///    u2(x,y) = \int u1(\xi,\eta) \delta(x-\xi) \delta(y-\eta) d\xi d\eta
    ///      \approx \sum u1(\xi,\eta) \delta(x-\xi) \delta(y-\eta) * dx^2
    /// The flux returned is the corresponding flux through cell edges,
    /// or u2 * dx
    Flux toFlux(const BoundaryVector& u) const;

    /// \brief Interpolate grid data to boundary.
    /// In particular, if q denotes fluxes in the 2D domain, compute
    /// corresponding velocities u2 = q / dx, and define velocities u1 at
    /// boundary points by interpolation:
    ///    u1(x,y) = \int u2(\xi,\eta) \delta(x-\xi) \delta(y-\eta) d\xi d\eta
    ///      \approx \sum u2(\xi,\eta) \delta(x-\xi) \delta(y-\eta) * dx^2
    BoundaryVector toBoundary(const Flux& u) const;

private:
    const Grid& _grid;
    const Geometry& _geometry;

    // Associations between BoundaryVector points and nearby Flux values
    struct Association {
        BoundaryVector::index boundaryIndex;
        Flux::index fluxIndex;
        double weight;
    };
    
    vector<Association> _neighbors;
};

} // namespace ibpm

#endif /* _REGULARIZER_H_ */
