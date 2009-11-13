#ifndef _OUTPUTTECPLOT_H_
#define _OUTPUTTECPLOT_H_

#include "Output.h"
#include "Scalar.h"
#include <string>
#include <vector>
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

// Small class for storing a list of pointers to variables, and their names
class VarList {
public:
    void addVariable( const Scalar* var, string name ) {
        _vars.push_back( var );
        _names.push_back( name );
    }
    
    int getNumVars() const {
        return _vars.size();
    }
    
    string getName(int i) const {
        return _names[i];
    }
    
    const Scalar* getVariable(int i) const {
        return _vars[i];
    }
    
private:
    vector<const Scalar*> _vars;
    vector<string> _names;
};
    
class OutputTecplot : public Output {
public:
    /// \brief Constructor
    /// \param[in] filename Filename in the standard printf format (e.g. "file%06d.plt"), where timestep will be supplied
    /// \param[in] title Title in the standard printf format
    OutputTecplot( string filename, string title );
    
    // Opens a new file with the specified name, and writes an ASCII Tecplot file
    // containing the variables in the VarList passed in
    bool writeTecplotFileASCII( const char* filename, const char* title, const VarList& list );
    
    /// \brief Write the Tecplot file
    bool doOutput(const State& x);
    
    /// \brief Change the filename for the output file
    void setFilename( string filename );
    
private:
    string _filename;
    string _title;
};

} // namespace ibpm

#endif /* _OUTPUTTECPLOT_H_ */
