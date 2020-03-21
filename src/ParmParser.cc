// ParmParser.cc
//
// Description:
// Implementation of class for parsing command-line arguments
//
// Author(s):
// Clancy Rowley
//
// Date: 4 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "ParmParser.h"
#include <fstream>
#include <iomanip>

const int OPTION_WIDTH = 20;

ParmParser::ParmParser(int argc, char* argv[]) :
    _argc( argc ),
    _argv( argc ),
    _used( argc ) {
    _args = "";
    for (int i=1; i<argc; ++i) {
        _args += argv[i];
        _args += " ";
        _argv[i] = argv[i];
        _used[i] = false;
    }
    _argv[0] = argv[0];
    _argOut << _argv[0];
    // remove path from executable name, if specified on command line
    string exec = _argv[0];
    int slashPosition = exec.find_last_of( '/' );
    if ( slashPosition >= 0 ) {
        exec = string( exec, slashPosition + 1, exec.length() );
    }
    _usageMessage << "USAGE: " << exec << " [options]" << endl
        << "where [options] are as follows (defaults shown in [ ]):" << endl;
}

void ParmParser::appendUsageMessage(
    string parm,
    string type,
    string description,
    string defaultVal
    ) {
    _usageMessage.setf(ios_base::left);
    _usageMessage << "  " << setw(OPTION_WIDTH) << "-" + parm + " " + type
        << description << " [" << defaultVal << "]" << endl;
}

void ParmParser::appendUsageMessage(
    string flag,
    string description
    ) {
    _usageMessage.setf(ios_base::left);
    _usageMessage << "  " << setw(OPTION_WIDTH) << "-" + flag
        << description << endl;
}

bool ParmParser::getFlag( string flag, string description ) {
    appendUsageMessage( flag, description );
    flag = "-" + flag;
    
    // Make an input stream with the whole argument list
    istringstream in;
    in.str(_args);
    
    // Loop through all the command line arguments
    int i=0;
    string s;
    while ( in >> s ) {
        ++i;
        // If the specified parameter is found
        if ( s == flag ) {
            _used[i] = true;
            _argOut << " " << flag;
            return true;
        }
    }
    // If not found, return false
    return false;
}

// Generic function to search for the given entry parm and return its
// argument or a default value
template<class T> T ParmParser::getParm(
    string parm,
    string description,
    T defaultVal
    ) {
    parm = "-" + parm;
    
    // Make an input stream with the whole argument list
    istringstream in;
    in.str(_args);
    
    // Loop through all the command line arguments
    int i=0;
    string s;
    while ( in >> s ) {
        ++i;
        // If the specified parameter is found
        if ( s == parm ) {
            // Try to get the next argument in the list as an integer
            T val;
            in >> val;
            // If successful
            if ( ! in.fail() ) {
                _used[i] = true;
                _used[i+1] = true;
                // Append the argument and value to the local list
                _argOut << " " << parm << " " << val;
                // return the value found
                return val;
            }
            else {
                // Otherwise, print an error message and return default value
                cerr << "Warning: cannot parse argument " << i 
                    << ": " << parm << endl;
                return defaultVal;
            }
        }
    }
    // If not found, return the default value
    return defaultVal;        
}

#define PP_APPEND_USAGE(a) \
    ostringstream defaultStr; \
    defaultStr << defaultVal; \
    appendUsageMessage( parm, a, description, defaultStr.str() );

// Search the parameter list for the given entry parm, returning
// defaultVal if not specified or invalid
int ParmParser::getInt( string parm, string description, int defaultVal ) {
    PP_APPEND_USAGE("<int>");
    return getParm( parm, description, defaultVal );
}

// Make a template function from getInt?
double ParmParser::getDouble( string parm, string description, double defaultVal ) {
    PP_APPEND_USAGE("<real>");
    return getParm( parm, description, defaultVal );
}

string ParmParser::getString( string parm, string description, string defaultVal ) {
    PP_APPEND_USAGE("<string>");
    return getParm( parm, description, defaultVal );
}

bool ParmParser::getBool( string parm, string description, bool defaultVal ) {
    PP_APPEND_USAGE("[0 or 1]");
    return getParm( parm, description, defaultVal );
}

#undef PP_APPEND_USAGE


bool ParmParser::inputIsValid() {
    // Loop through all the parameters, checking if have been used
    bool allused = true;
    for (int i=1; i<_argc; ++i) {
        allused = allused && _used[i];
    }

    if (allused) {
        return true;
    }
    else {
        cerr << "Warning: the following parameters were not used:"
            << endl << "  ";
        for (int i=1; i<_argc; ++i) {
            if (! _used[i]) {
                cerr << _argv[i] << " ";
            }
        }
        cerr << endl;
        return false;
    }
}


void ParmParser::printUsage( ostream& out ) {
    out << _usageMessage.str();
}

string ParmParser::getParameters() {
    return _argOut.str();
}

// Save argument list to a file
void ParmParser::saveParameters( string fname ) {
    fstream out( fname.c_str(), ios_base::out );
    out << _argOut.str() << endl;
}
