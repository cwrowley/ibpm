#ifndef _OUTPUTPROBES_H_
#define _OUTPUTPROBES_H_

#include "Output.h"
#include "Grid.h"
#include "Array.h"
#include <vector>
#include <string>
using namespace std;

namespace ibpm {

/*!
\file OutputProbes.h
\class OutputProbes

\brief Write velocities, fluxes, and vorticity, at given probe locations, to files. 
 
Each probe has a corresponding output file. 

All probes are supposed to be located at the interior nodes
 at the finest grid level (level 0).

Probes are labelled as Probe 1, 2, ... .

Probe information (probe #, position) is stored in a separate file.

\author Clancy Rowley
\author Zhanhua Ma
\date 11 May 2009

*/

class OutputProbes : public Output {
public:
	/// \brief Constructor
    /// \param[in] filename, to which probe data will be written.
	OutputProbes(string filename, Grid& grid);

    /// Open the file with name (filename + "info"), 
	/// and write probe information (probe #, probe position).
	/// Also, for each probe, open a file with name (filename + "probe #").
    /// If a file with the same name is already present, it is overwritten.
    /// Returns true if successful.
    bool init();

    /// \brief Close all the files.
    /// Returns true if successful
    bool cleanup(); 
	
    /// Write velocities u, v, fluxes q.x, q,y and vorticity omega for each probe, 
	/// to the correpsonding file with name (filename + probe#).
	bool doOutput( const State& x );
	
	/// \brief Add a probe by specifying its gridpoint indices 
	void addProbeByIndex( int i, int j );
	
	/// \brief Add a probe by specifying its absolute coordinates 
	void addProbeByPosition( double xcord, double ycord );
	// TODO: Write up this member function.	
	
	/// \brief Add a probe by specifying its gridpoint indices
	void addProbe( int i, int j );
	
	/// \brief Add a probe by specifying its absolute coordinates
	void addProbe( double xcord, double ycord ); 
	
	/// Print out probe locations (by grid indices), for debugging 
	void print();
	
	/// Return the number of probes
	inline int getNumProbes(){ return _numProbes; }
	
	/// Return the gridpoint index i of the corresponding probe
	inline int getProbeIndexX( int index ){
		assert( index <= _numProbes && index >= 1 ); 
		int i = _probePositions( index - 1, 0 );
		return i;
	}
	
	/// Return the gridpoint index j of the corresponding probe
	inline int getProbeIndexY( int index ){
		assert( index <= _numProbes ); 
		assert( index >= 1 );
		int j = _probePositions( index - 1, 1 );
		return j;
	}
	
	/// Return the gridpoint x coordinate of the corresponding probe
	inline double getProbeCoordX( int index ){
		assert( index <= _numProbes && index >= 1 ) ; 
		double xcoord = _grid.getXEdge( _lev, getProbeIndexX( index )); 
		return xcoord;
	}
	
	/// Return the gridpoint y coordinate of the corresponding probe
	inline double getProbeCoordY( int index ){
		assert( index <= _numProbes && index >= 1 ); 
		double ycoord = _grid.getYEdge( _lev, getProbeIndexY( index )); 
		return ycoord;
	}
	
private:
    string _filename;
	Grid _grid;
	FILE* _fp;
	int _numProbes;
	bool _flagInitialization;  //check if has tried initialization 
	int _lev;
	int _dimen; 
	Array::Array2<int> _probePositions; // probe positions, represented by grid indices
	vector<FILE*> _fprobes; 
};

} // namespace ibpm

#endif /* _OUTPUTPROBES_H_ */
