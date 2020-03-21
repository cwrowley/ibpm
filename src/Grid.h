#ifndef _GRID_H_
#define _GRID_H_

#include <assert.h>
#include <math.h>
//JON
#include <iostream>

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

    /// Constructor: set all standard grid parameters
    Grid(
         int nx,
         int ny,
         int ngrid,
         double length,
         double xOffset,
         double yOffset
    );
	
	/// Constructor: set all grid parameters, including x- and y- shifts
    Grid(
         int nx,
         int ny,
         int ngrid,
         double length,
         double xOffset,
         double yOffset,
		 double xShift,
		 double yShift
		 );
	
    /// Default constructor: set all parameters to zero
    Grid();

    /// Set all grid parameters
    void resize(
        int nx,
        int ny,
        int ngrid,
        double length,
        double xOffset,
        double yOffset
    );
	
	/// Set all grid parameters, including x- and y- shifts
    void resize(
		int nx,
		int ny,
		int ngrid,
		double length,
		double xOffset,
		double yOffset,
		double xShift,
		double yShift
	);
    /// Return number of cells in x-direction
    inline int Nx() const { return _nx; }

    /// Return number of coarse cells outside each fine domain, in x-direction
	// This counts the number of cells to the left of the fine domain
	// This distinction is key when the grid is shifted
	// To round to the nearest integer we add 0.5 and use the floor function
	// However, due to asserts enforced in Grid.cc, this is nominally an integer already
    inline int NxExt() const { 
		return (int) floor( _nx / 4 * ( 1 - _xShift ) + 0.5 ); 
	}
    
    /// Return number of cells in y-direction
    inline int Ny() const { return _ny; }
    
    /// Return number of coarse cells outside each fine domain, in y-direction
	// This counts the number of cells below the fine domain
	// This distinction is key when the grid is shifted
	// To round to the nearest integer we add 0.5 and use the floor function
	// However, due to asserts enforced in Grid.cc, this is nominally an integer already
    inline int NyExt() const { return (int) floor( _ny / 4 * ( 1 - _yShift ) + 0.5 ); }
    
    /// Given indices (i,j) on a coarse grid, return corresponding indices
    /// (ii,jj) on the fine grid
    inline void c2f(int i, int j, int& ii, int& jj) const {
        ii = (i - NxExt() ) * 2;
        jj = (j - NyExt() ) * 2;
    }

    /// Given indices (ii,jj) on a fine grid, return indices (i,j) of the
    /// corresponding point on the coarse grid, or of the nearest point
    /// below and to the left
    inline void f2c(int ii, int jj, int& i, int& j) const {
        i = ii/2 + NxExt();
        j = jj/2 + NyExt();
    }
    
    /// Return number of grid levels for multi-domain solution
    inline int Ngrid() const { return _ngrid; }
    
    /// Return grid spacing on finest level (same in x- and y-directions)
    inline double Dx() const { return _dx; }
    
    /// Return grid spacing at specified grid level
    inline double Dx(int lev) const {
        assert( lev >= 0 && lev < _ngrid );
        return _dx * (1 << lev);
    }

    /// Return the x-coordinate of the center of cell i  (i in 0..m-1)
    double getXCenter(int lev, int i) const;
    
    /// Return the y-coordinate of the center of cell j  (j in 0..n-1)
    double getYCenter(int lev, int j) const;

    /// Return the x-coordinate of the left edge of cell i  (i in 0..m)
    double getXEdge(int lev, int i) const;
    
    /// Return the y-coordinate of the bottom edge of cell j  (j in 0..n)
    double getYEdge(int lev, int j) const;
	
    /// Return the grid index i corresponding to the given x-coordinate 
	/// Currently, only works for the finest grid level. 
    int getXGridIndex( double x ) const;
	
    /// Return the grid index j corresponding to the given x-coordinate 
	/// Currently, only works for the finest grid level. 
    int getYGridIndex( double y ) const;
	
    /// Set shift parameter in x
    void setXShift(double xShift);

    /// Get the current x-shift parameter
    double getXShift() const;
	
	/// Set shift parameter in y
    void setYShift(double yShift);
	
    /// Get the current y-shift parameter
    double getYShift() const;
    
    /// Compare two grids
    bool isEqualTo( const Grid& grid2 ) const;

private:
    double getXOffset( int lev ) const;
    double getYOffset( int lev ) const;
    int _nx;
    int _ny;
    int _ngrid;
    double _dx;
    double _xOffset;
    double _yOffset;
    double _xShift;
	double _yShift;
};

} // namespace

#endif /* _GRID_H_ */
