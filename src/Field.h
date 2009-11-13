#ifndef _FIELD_H_
#define _FIELD_H_

#include "Grid.h"

namespace ibpm {

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

class Field {
public:
    Field();
    Field( const Grid& grid );
    Field( const Field& field );
    virtual ~Field();
    
    inline int Nx() const { return _grid.Nx(); }
    inline int NxExt() const { return _grid.NxExt(); }
    inline int Ny() const { return _grid.Ny(); }
    inline int NyExt() const { return _grid.NyExt(); }
    inline int Ngrid() const { return _grid.Ngrid(); }
    inline double Dx() const { return _grid.Dx(); }
    inline double Dx(int lev) const { return _grid.Dx(lev); }
    inline double getXShift() const { return _grid.getXShift(); }
    inline double getYShift() const { return _grid.getYShift(); }
    inline const Grid& getGrid() const { return _grid; }
    inline void setGrid( const Grid& grid ) { _grid = grid; }

    /// Return the x-coordinate of the center of cell i  (i in 0..m-1)
    inline double getXCenter(int lev, int i) const {
        return _grid.getXCenter( lev, i );
    }
    
    /// Return the y-coordinate of the center of cell j  (j in 0..n-1)
    inline double getYCenter(int lev, int j) const {
        return _grid.getYCenter( lev, j );
    }

    /// Return the x-coordinate of the left edge of cell i  (i in 0..m)
    inline double getXEdge(int lev, int i) const {
        return _grid.getXEdge( lev, i );
    }
    
    /// Return the y-coordinate of the bottom edge of cell j  (j in 0..n)
    inline double getYEdge(int lev, int j) const {
        return _grid.getYEdge( lev, j );
    }
    
private:
    Grid _grid;
};

} // namespace ibpm


#endif /* _FIELD_H_ */
