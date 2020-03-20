// utils.cc
//
// Description:
// Utilities for working with strings
//
// Author(s):
// Clancy Rowley
//
// Date: 1 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include <string>
#include <ctype.h>

using std::string;

void EatWhitespace( string& s ) {
    while ( isspace( s[0] ) ) {
        s.erase(0,1);
    }
}

void MakeLowercase( string& s ) {
    for (unsigned int i=0; i<s.length(); ++i) {
        s[i] = tolower( s[i] );
    }
}

void AddSlashToPath( string& s) {
    int len = s.length();
    if ( len > 0 && s[len-1] != '/' ) {
        s += "/";
    }
}
