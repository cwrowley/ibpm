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
    /// \param[in] formatString Filename in the standard printf format (e.g. "file%06d.bin")
    OutputRestart(string formatString);

    /// \brief Write the restart file
    bool doOutput(const State& x);
    
private:
    string _formatString;
};

} // namespace ibpm

#endif /* _OUTPUTRESTART_H_ */
