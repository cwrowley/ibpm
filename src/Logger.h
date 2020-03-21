#ifndef _LOGGER_H_
#define _LOGGER_H_

#include "BaseFlow.h"
#include "State.h"
#include <vector>

using std::vector;

namespace ibpm {

class Output;

/*!
    \file Logger.h
    \class Logger
    
    \brief Maintain a list of output routines, and call them when specified.
    
    \author Clancy Rowley
    \date 21 Aug 2008

    \author $LastChangedBy$
    \date $LastChangedDate$
    \version $Revision$
*/

class Logger {
public:
    Logger();
    
    /// \brief Add the specified Output to the list of output routines.
    /// The caller is responsible for allocating and deallocating memory for
    /// Output.
    /// \param[in] output is a pointer to the output routine to add.
    /// \param[in] numSkip specifies how often to call the given routine
    ///             (number of timesteps).
    void addOutput(Output* output, int numSkip);

    /// \brief Call all output routines needed at the current timestep.
    bool doOutput(const State& x);
    
    /// \brief Call all output routines needed at the current timestep, making use of baseflow.
    bool doOutput(const BaseFlow& q, const State& x);
    
    /// \brief Initialize all of the output routines.
    bool init();

    /// \brief Clean up all of the output routines.
    bool cleanup();

private:
    struct Entry {
        Output* output;
        int numSkip;
        inline bool shouldBeCalled(const State& x) {
            return (x.timestep % numSkip == 0);
        }
    };
    vector<Entry> _outputs;
    bool _hasBeenInitialized;
};

#undef LOOP_OVER_ALL_OUTPUTS

} // namespace ibpm

#endif /* _LOGGER_H_ */
