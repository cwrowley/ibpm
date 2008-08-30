#ifndef _OUTPUTTECPLOT_H_
#define _OUTPUTTECPLOT_H_

#include "Output.h"
#include <string>
using std::string;

namespace ibpm {

/*!
\file OutputTecplot.h
\class OutputTecplot

\brief Output routines for writing ASCII Tecplot files

\author Clancy Rowley
\date 22 Aug 2008

\author $LastChangedBy$
\date $LastChangedDate$
\version $Revision$
*/

class OutputTecplot : public Output {
public:
    /// \brief Constructor
    /// \param[in] filename Filename in the standard printf format (e.g. "file%06d.plt"), where timestep will be supplied
    /// \param[in] title Title in the standard printf format
    OutputTecplot( string filename, string title );

    /// \brief Write the Tecplot file
    bool doOutput(const State& x);
    
private:
    string _filename;
    string _title;
};

} // namespace ibpm

#endif /* _OUTPUTTECPLOT_H_ */
