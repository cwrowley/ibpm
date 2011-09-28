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

\author $LastChangedBy$
\date $LastChangedDate$
\version $Revision$
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
    
    /// \brief Write data to the force file from a state object, angle of attack, and freestream velocity.
    bool doOutput(const double alpha, const double mag, const State& x);

    /// \brief Write data to the force file, making use of baseflow.
    bool doOutput(const BaseFlow& q, const State& x);
    
    /// \brief Write data to the force file.
    bool doOutput(const State& x);
    
    
private:   
    string _filename;
    FILE* _fp;
};

} // namespace ibpm

#endif /* _OUTPUTFORCE_H_ */
