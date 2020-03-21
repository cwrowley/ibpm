#ifndef _STATE_H_
#define _STATE_H_

#include "Grid.h"
#include "Flux.h"
#include "Scalar.h"
#include "BoundaryVector.h"
#include <string>

using std::string;

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
    /// Default constructor: do not allocate memory
    State();
    	  		
    State( const Grid& grid, int numPoints );

    /// \brief Instantiate a state by reading data from the specified file
    State( string filename );

    ~State();
    
    /// \brief Allocate memory, with the specified Grid and number of
    /// boundary points
    void resize( const Grid& grid, int numPoints );

    /// \brief Save the state to a file (e.g. as a restart file)
    /// Return true if successful
    /*  WARNING:  At this point, the xshift and yshift parameters are not saved
     and are not checked for compatibility when loading.  Caution should be
     taken when working with shifted grids.  This approach was taken to pre-
     serve backwards compatibility with previously saved binary files.  In 
     the future perhaps using HDF5 would prevent such problems.
     */
    bool save(std::string filename) const;

    /// \brief Load the state from a file (e.g. as a restart file)
    /// Return true if successful
    bool load(const std::string& filename);

    /// \brief Routine for computing X & Y forces
    void computeNetForce( double& xforce, double& yforce ) const;
			    
    // public data
    Flux q;
    Scalar omega;
    BoundaryVector f;
    int timestep;
    double time;
};

} // namespace ibpm

#endif /* _STATE_H_ */
