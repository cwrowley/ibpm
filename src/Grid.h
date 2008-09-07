#ifndef _GRID_H_
#define _GRID_H_

namespace ibpm {

/*!
    \file Grid.h
    \class Grid

    \brief Define parameters associated with a uniform staggered grid

    Uses a staggered grid, suitable for a finite-volume method, in which
    scalars (such as streamfunction and vorticity) are defined at cell nodes
    and fluxes are defined at cell edges.  (Note that pressure is normally
    defined at cell centers, but pressure is not needed in this simulation.)

    \image html grid.pdf "Layout of staggered grid"
    \image latex grid.pdf "Layout of staggered grid"

    \author Clancy Rowley
    \author $LastChangedBy$
    \date  3 Jul 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class Grid {
public:
    /// Constructor: set all grid parameters
    Grid( int nx, int ny, double length, double xOffset, double yOffset );

    /// Default constructor: set all parameters to zero
    Grid();

    /// Set all grid parameters
    void resize(
        int nx,
        int ny,
        double length,
        double xOffset,
        double yOffset
    );
    
    /// Return number of cells in x-direction
    inline int getNx() const { return _nx; }
    
    /// Return number of cells in y-direction
    inline int getNy() const { return _ny; }
    
    /// Return grid spacing (same in x- and y-directions)
    inline double getDx() const { return _dx; }

    /// Return the x-coordinate of the center of cell i  (i in 0..m-1)
    double getXCenter(int i) const;
    
    /// Return the y-coordinate of the center of cell j  (j in 0..n-1)
    double getYCenter(int j) const;

    /// Return the x-coordinate of the left edge of cell i  (i in 0..m)
    double getXEdge(int i) const;
    
    /// Return the y-coordinate of the bottom edge of cell j  (j in 0..n)
    double getYEdge(int j) const;
    
    /// Return the area of the domain
    double getArea() const;
    

private:
    int _nx;
    int _ny;
    double _dx;
    double _xOffset;
    double _yOffset;
};

} // namespace

#endif /* _GRID_H_ */
