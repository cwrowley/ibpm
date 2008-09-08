#ifndef _FIELD_H_
#define _FIELD_H_

#include "Grid.h"

/*!
    \file Field.h
    \class Field

    \brief Abstract base class for a field defined on a grid.
    
    For instance, scalar fields (Scalar) and vector fields described by
    Flux objects are subclasses of a Field.

    \author Clancy Rowley
    \author $LastChangedBy$
    \date  7 Sep 2008
    \date $LastChangedDate$
    \version $Revision$
*/

namespace ibpm {

class Field {
public:
    Field();
    Field( const Grid& grid );
    virtual ~Field();
    
    inline int getNx() const { return _grid.getNx(); }
    inline int getNy() const { return _grid.getNy(); }
    inline double getDx() const { return _grid.getDx(); }
    inline const Grid& getGrid() const { return _grid; }
    inline void setGrid( const Grid& grid ) { _grid = grid; }

    /// Return the x-coordinate of the center of cell i  (i in 0..m-1)
    inline double getXCenter(int i) const { return _grid.getXCenter( i ); }
    
    /// Return the y-coordinate of the center of cell j  (j in 0..n-1)
    inline double getYCenter(int j) const { return _grid.getYCenter( j ); }

    /// Return the x-coordinate of the left edge of cell i  (i in 0..m)
    inline double getXEdge(int i) const { return _grid.getXEdge( i ); }
    
    /// Return the y-coordinate of the bottom edge of cell j  (j in 0..n)
    inline double getYEdge(int j) const { return _grid.getYEdge( j ); }
    
    /// Return the area of the domain
    inline double getArea() const { return _grid.getArea(); }

private:
    Grid _grid;
};

} // namespace ibpm


#endif /* _FIELD_H_ */
