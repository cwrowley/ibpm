#ifndef _DIRECTION_H_
#define _DIRECTION_H_

namespace ibpm {

/*!
    \file Direction.h

    \brief Type for specifying directions of vectors (X or Y)

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 27 Jul 2008
    \date $LastChangedDate$
    \version $Revision$
*/
enum Direction {X, Y, XY};

/// Postfix operator dir++
inline Direction operator++(Direction& dir, int) {
    return dir = (Direction)(dir + 1);
}

/// Prefix operator ++dir
inline Direction operator++(Direction& dir) {
    return dir = (Direction)(dir + 1);
}

} // namespace ibpm

#endif /* _DIRECTION_H_ */
