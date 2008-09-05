#ifndef _PARMPARSER_H_
#define _PARMPARSER_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
using namespace std;

/*!
    \file ParmParser.h
    \class ParmParser

    \brief Class for parsing command-line arguments

    \author Clancy Rowley
    \author $LastChangedBy$
    \date  4 Sep 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class ParmParser {
public:
    /// Constructor, taking command line arguments as input
    ParmParser(int argc, char* argv[]);

    /// Search the parameter list for the given flag, that does not take
    /// an argument.  Returns true if the flag is present
    bool getFlag( string flag, string description );
    
    /// Search the parameter list for the given entry parm and a single
    /// argument of the indicated type, returning defaultVal if not specified
    /// If argument is invalid, print a warning message and return defaultVal
    int getInt( string parm, string description, int defaultVal );
    double getDouble( string parm, string description, double defaultVal );
    string getString( string parm, string description, string defaultVal );
    bool getBool( string parm, string description, bool defaultVal );

    /// \brief Check if any input parameters were invalid, or unused
    /// Returns true if everything is valid and all parameters used
    /// If any were not used, prints them to standard error and returns false
    bool inputIsValid();

    /// \brief Print a usage message
    void printUsage( ostream& out );
    
    /// \brief Return argument list in a string
    string getParameters();
    
    /// \brief Save argument list to a file
    void saveParameters( string fname );

private:
    void appendUsageMessage( string flag, string description );
    void appendUsageMessage(
        string parm,
        string type,
        string description,
        string defVal
    );
    template<class T> T getParm(
        string parm,
        string description,
        T defaultVal
    );
    
    // data
    int _argc;
    vector<string> _argv;
    vector<bool> _used;
    string _args;
    ostringstream _argOut;
    ostringstream _usageMessage;
};

#endif /* _PARMPARSER_H_ */
