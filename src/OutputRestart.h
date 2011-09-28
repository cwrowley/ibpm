#ifndef _OUTPUTRESTART_H_
#define _OUTPUTRESTART_H_

#include "Output.h"
#include <string>
using std::string;

namespace ibpm {

/*!
\file OutputRestart.h
\class OutputRestart

\brief Output routine for writing a restart file

\author Clancy Rowley
\date 21 Aug 2008

\author $LastChangedBy$
\date $LastChangedDate$
\version $Revision$
*/

class OutputRestart : public Output {
public:
    /// \brief Constructor
    /// \param[in] formatString Filename in the standard printf format
    /// (e.g. "file%06d.bin"), where timestep will be substituted for %d
    OutputRestart(string formatString);

    /// \brief Write the restart file
    bool doOutput(const State& x);
    
    /// \brief Write the restart file, making use of the baseflow too
    bool doOutput(const BaseFlow& q, const State& x);
    
    /// \brief Change the filename for the output file
    void setFilename( string formatString );
    
private:
    string _formatString;
};

} // namespace ibpm

#endif /* _OUTPUTRESTART_H_ */
