// Logger.cc
//
// Description:
// Implementation of Logger class for calling output routines
//
// Author(s):
// Clancy Rowley
//
// Date: 22 Aug 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include <assert.h>
#include <vector>
#include "Output.h"
#include "BaseFlow.h"
#include "State.h"
#include "Logger.h"


namespace ibpm {

#define LOOP_OVER_ALL_OUTPUTS( a )                                          \
    vector<Entry>::iterator entry;                                          \
    bool successful = true;                                                 \
    for (entry = _outputs.begin(); entry != _outputs.end(); ++entry) {      \
        bool result = (entry->output)->a;                                   \
        successful = successful && result;                                  \
    }                                                                       \
    return successful;


Logger::Logger() {
    _hasBeenInitialized = false;
}

void Logger::addOutput(Output* output, int numSkip) {
    Entry entry;
    entry.output = output;
    entry.numSkip = numSkip;
    _outputs.push_back(entry);
}

/// \brief Call all output routines needed at the current timestep,.
bool Logger::doOutput(const State& x) {
    assert( _hasBeenInitialized );
    vector<Entry>::iterator entry;
    bool successful = true;
    for (entry = _outputs.begin(); entry != _outputs.end(); ++entry) {
        if ( entry->shouldBeCalled( x ) ) {
            bool result = (entry->output)->doOutput( x );
            successful = successful && result;
        }
    }
    return successful;
}
    
/// \brief Call all output routines needed at the current timestep, using baseflow.
bool Logger::doOutput(const BaseFlow& q, const State& x) {
	assert( _hasBeenInitialized );
    vector<Entry>::iterator entry;
    bool successful = true;
    for (entry = _outputs.begin(); entry != _outputs.end(); ++entry) {
        if ( entry->shouldBeCalled( x ) ) {
            bool result = (entry->output)->doOutput( q , x );
            successful = successful && result;
        }
    }
    return successful;
}

// Initialize all of the output routines.
bool Logger::init() {
    _hasBeenInitialized = true;
    LOOP_OVER_ALL_OUTPUTS( init() );
}

/// Clean up all of the output routines.
bool Logger::cleanup() {
    assert( _hasBeenInitialized );
    LOOP_OVER_ALL_OUTPUTS( cleanup() );
}

#undef LOOP_OVER_ALL_OUTPUTS

} // namespace
