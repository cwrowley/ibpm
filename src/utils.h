#ifndef _UTILS_H_
#define _UTILS_H_

/*!
    \file utils.h

    \brief Utilities for working with strings

    \author Clancy Rowley
    \author $LastChangedBy$
    \date  2 Sep 2008
    \date $LastChangedDate$
    \version $Revision$
*/

/// \brief Remove whitespace at the beginning of the string s
void EatWhitespace( string& s );

/// \brief Convert the string s to lower case
void MakeLowercase( string& s );

/// \brief Add a slash to the end of the string s, if not already present
void AddSlashToPath( string& s );

#endif /* _UTILS_H_ */
