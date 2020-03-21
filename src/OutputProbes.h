#ifndef _OUTPUTPROBES_H_
#define _OUTPUTPROBES_H_

#include "Output.h"
#include "Grid.h"
#include "Array.h"
#include <vector>
#include <string>

using std::string;
using std::vector;

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
    /// For instance, if filename = "out/probe%02.dat"
    /// then the following files will be created:
    ///   out/probe00.dat   (description of probe locations)
    ///   out/probe01.dat   (first probe)
    ///   out/probe02.dat   (second probe)
    ///   ...
    OutputProbes(string filename, Grid& grid);

    /// Write a file with description of probe locations
    ///   (this file has index 0)
    /// Also, open a file for each probe
    /// If a file with the same name is already present, it is overwritten.
    /// Returns true if successful.
    bool init();

    /// \brief Close all the files.
    /// Returns true if successful
    bool cleanup(); 
    
    /// Write velocities u, v, fluxes q.x, q,y and vorticity omega for each probe, 
    /// to the correpsonding file with name (filename + probe#).
    bool doOutput( const State& x );
    
    /// Write velocities u, v, fluxes q.x, q,y and vorticity omega for each probe, 
    /// to the correpsonding file with name (filename + probe#), making use of the baseflow too.
    bool doOutput( const BaseFlow& q, const State& x );
    
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
    inline int getNumProbes(){ return _probes.size(); }
    
    /// Return the gridpoint index i of the corresponding probe
    inline int getProbeIndexX( unsigned int index ){
        assert( index <= _probes.size() && index >= 1 ); 
        return _probes[index-1].i;
    }
    
    /// Return the gridpoint index j of the corresponding probe
    inline int getProbeIndexY( unsigned int index ){
        assert( index <= _probes.size() ); 
        assert( index >= 1 );
        return _probes[index-1].j;
    }
    
    /// Return the gridpoint x coordinate of the corresponding probe
    inline double getProbeCoordX( unsigned int index ){
        assert( index <= _probes.size() && index >= 1 );
        return _grid.getXEdge( _lev, _probes[index-1].i );
    }
    
    /// Return the gridpoint y coordinate of the corresponding probe
    inline double getProbeCoordY( unsigned int index ){
        assert( index <= _probes.size() && index >= 1 );
        return _grid.getYEdge( _lev, _probes[index-1].j );
    }
    
private:
    class Probe {
    public:
        Probe(int ii, int jj) :
            i(ii),
            j(jj),
            fp(NULL)
        {}
        int i,j;
        FILE *fp;
    };

    bool writeSummaryFile(void);

    // Private data
    string _filename;
    Grid _grid;
    bool _hasBeenInitialized;
    vector<Probe> _probes;
    static const int _lev;
    static const int _dimen;
};

} // namespace ibpm

#endif /* _OUTPUTPROBES_H_ */
