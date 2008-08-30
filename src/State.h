#ifndef _STATE_H_
#define _STATE_H_

#include "Grid.h"
#include "Geometry.h"
#include "Flux.h"
#include "Scalar.h"
#include "BoundaryVector.h"

namespace ibpm {

/*!
    \file State.h
    \class State

    \brief Structure for grouping state variables

    \author Clancy Rowley
    \author $LastChangedBy$
    \date  7 Jul 2008
    \date $LastChangedDate$
    \version $Revision$
*/

struct State {
    State(const Grid& grid, const Geometry& geom) :
        q(grid),
        gamma(grid),
        f( geom.getNumPoints() ),
        timestep(0),
        time(0.)
    {}

    ~State() {}
    
    void loadRestartFile(string filename) {}
    
    Flux q;
    Scalar gamma;
    BoundaryVector f;
    int timestep;
    double time;
};

} // namespace ibpm

#endif /* _STATE_H_ */
