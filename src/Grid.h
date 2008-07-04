#ifndef _GRID_H_
#define _GRID_H_

/*!
	\file Grid.h
	\class Grid

	\brief Define parameters associated with a uniform staggered grid

	Uses a staggered grid, suitable for a finite-volume method, in which
	scalars are defined at cell centers and fluxes are defined at cell edges.
	
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
	Grid(int nx, int ny, double length, double xOffset, double yOffset) :
		_nx(nx),
		_ny(ny),
		_xOffset(xOffset),
		_yOffset(yOffset)
	{
		_dx = length / nx;
	}

	/// Return number of cells in x-direction
	inline int getNx() { return _nx; }
	
	/// Return number of cells in y-direction
	inline int getNy() { return _ny; }
	
	/// Return grid spacing (same in x- and y-directions)
	inline double getDx() { return _dx; }

	/// Return the x-coordinate of the center of cell i  (i in 0..m-1)
	double getXCenter(int i);
	
	/// Return the y-coordinate of the center of cell j  (j in 0..n-1)
	double getYCenter(int j);

	/// Return the x-coordinate of the left edge of cell i  (i in 0..m)
	double getXEdge(int i);
	
	/// Return the y-coordinate of the bottom edge of cell j  (j in 0..n)
	double getYEdge(int j);


private:
	int _nx;
	int _ny;
	double _dx;
	double _xOffset;
	double _yOffset;
};

#endif /* _GRID_H_ */
