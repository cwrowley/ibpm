#ifndef _STATE_H_
#define _STATE_H_

#include "Grid.h"
#include "Geometry.h"
#include "Flux.h"
#include "Scalar.h"
#include "BoundaryVector.h"
#include <string>

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

class State {
public:
    State(const Grid& grid, const Geometry& geom);

    ~State();
    
    /// \brief Save the state to a file (e.g. as a restart file)
    /// Return true if successful
    bool save(std::string filename) const;

    /// \brief Load the state from a file (e.g. as a restart file)
    /// Return true if successful
    bool load(const std::string& filename);

    // public data
    Flux q;
    Scalar gamma;
    BoundaryVector f;
    int timestep;
    double time;
};

} // namespace ibpm

#endif /* _STATE_H_ */
