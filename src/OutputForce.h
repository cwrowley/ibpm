#ifndef _OUTPUTFORCE_H_
#define _OUTPUTFORCE_H_

#include "Output.h"
#include <string>
using std::string;

namespace ibpm {

/*!
\file OutputForce.h
\class OutputForce

\brief Output routine for writing a list of force coefficients

\author Steve Brunton
\date 21 Aug 2008

\author $LastChangedBy: cwrowley $
\date $LastChangedDate: 2008-09-01 21:43:39 -0400 (Mon, 01 Sep 2008) $
\version $Revision: 117 $
*/

class OutputForce : public Output {
public:
    /// \brief Constructor
    /// \param[in] filename Filename in the standard printf format
    /// (e.g. "file%06d.bin"), where timestep will be substituted for %d
    OutputForce(string filename);

    /// \brief Write the force file
    bool doOutput(const State& x);
    
private:
    string _filename;
    char _fname[256];
};

} // namespace ibpm

#endif /* _OUTPUTFORCE_H_ */
