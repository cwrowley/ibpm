#ifndef _OUTPUTFORCE_H_
#define _OUTPUTFORCE_H_

#include "Output.h"
#include <stdio.h>
#include <string>
using std::string;

namespace ibpm {

/*!
\file OutputForce.h
\class OutputForce

\brief Output routine for writing a list of force coefficients.

\author Steve Brunton
\author Clancy Rowley
\date 21 Aug 2008

\author $LastChangedBy: cwrowley $
\date $LastChangedDate: 2008-09-01 21:43:39 -0400 (Mon, 01 Sep 2008) $
\version $Revision: 117 $
*/

class OutputForce : public Output {
public:
    /// \brief Constructor
    /// \param[in] filename to which force data will be written.
    OutputForce(string filename);

    /// \brief Open the file for writing.
    /// If a file with the same name is already present, it is overwritten.
    /// Returns true if successful.
    bool init();

    /// \brief Close the file.
    /// Returns true if successful
    bool cleanup();

    /// \brief Write data to the force file.
    bool doOutput(const State& x);
    
private:
    string _filename;
    FILE* _fp;
};

} // namespace ibpm

#endif /* _OUTPUTFORCE_H_ */
