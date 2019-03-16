// bin2plt - utility to convert restart files to Tecplot format, for IBPM code
//
// Clancy Rowley
// 6 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include <string>
#include <iostream>
#include <fstream>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

string GetBasename( string s );

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <filenames>" << endl;        
        exit(1);
    }
    
    for (int i=1; i<argc; ++i) {
        State x;
        // Read in a restart file
        if ( ! x.load( argv[i] ) ) {
            cerr << "Error reading file " << argv[i] << endl;
        }
        else {
            // Write out the Tecplot file
            string outname = GetBasename( argv[i] );
            outname += ".plt";
            string title( argv[i] );
            OutputTecplot tecplot( outname, title, 1 );
            if ( ! tecplot.doOutput( x ) ) {
                cerr << "Error writing file " << outname << endl;
            }
        }
    }
    
    return 0;
}

// remove the leading part of the path, and the suffix
//   e.g. if s = "/path/to/myfile.bin"
//   returns "myfile"
string GetBasename( string s ) {
    int slashPosition = s.find_last_of( '/' );
    if ( slashPosition >= 0 ) {
        s = string( s, slashPosition + 1, s.length() );
    }
    int suffixPosition = s.find_last_of( '.' );
    if ( suffixPosition >= 0 ) {
        s = string( s, 0, suffixPosition );
    }
    return s;
}
