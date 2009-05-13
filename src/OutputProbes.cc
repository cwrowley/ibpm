// OutputProbes.cc
//
// Description:
// Implementation of output routine for writing a list of probes' states
//
// Author(s):
// Clancy Rowley
// Zhanhua Ma
//
// Date: 11 May 2008

#include "OutputProbes.h"
#include "Array.h"
#include "State.h"
#include "Output.h"
#include "VectorOperations.h"
#include <vector>
#include <stdio.h>
#include <string>
using namespace std;

namespace ibpm {

OutputProbes::OutputProbes(string filename, Grid& grid ) :
    _filename( filename ),
	_grid( grid ) {
		_numProbes = 0;
		_lev = 0; // All probes are located at the finest level of grid
		_dimen = 2; // Two dimensional domain for now.
		_flagInitialization = false; 
}

bool OutputProbes::init() {
	_flagInitialization = true;
    _fp = fopen( (_filename + "info").c_str(), "w" );
	if ( _fp == NULL ) return false;
	
	// write: probe#, probe grid indices i,  j,  probe coordinates x, y 
	assert( unsigned ( _numProbes * _dimen ) == _probePositions.Size() );
	for ( int i = 0 ; i < _numProbes ; i ++ ) {
		int gridi = getProbeIndexX( i + 1 );
		int gridj = getProbeIndexY( i + 1 );
		double xcord = getProbeCoordX( i + 1 );
		double ycord = getProbeCoordY( i + 1 );
		fprintf( _fp, "%2d %3d %d %.5e %.5e\n", 
				i + 1,  gridi, gridj, xcord, ycord ); 
	}
	fclose( _fp ); 
	
	// open files for each probe to output data
	FILE* pinit;
	_fprobes.assign( _numProbes - 1 ,  pinit ); 
	for ( int i = 0 ; i < _numProbes ; i ++ ) {
		char name[256];
		sprintf( name, (_filename + "%02d").c_str(), i + 1 );
		_fprobes[i] = fopen( name, "w" );
		if ( _fprobes[i] == NULL ) return false;
	}				
	    
    return true;
}

bool OutputProbes::cleanup() {
    bool status = true;
    if ( _fp != NULL ) {
        status = fclose( _fp );
    }
	for ( int i = 0 ; i < _numProbes ; i ++ ) {
		if ( _fprobes[i] != NULL ) {
		        status = fclose( _fprobes[i] );
		}
	}
    return status;
}

bool OutputProbes::doOutput( const State& state) {
	for ( int i = 0 ; i < _numProbes ; i ++ ) {
		if ( _fprobes[i] == NULL ) return false;
	}
	
	Scalar u(_grid);
    Scalar v(_grid);
	FluxToVelocity( state.q, u, v );

	// Write  u, v, qx, qy, omega, all at gridpoint/edge (i, j), for each probe,
	// in seperate files
	for ( int i = 0; i < _numProbes ; i ++ ) {
//		int gridi = _probePositions( i, 0 );
//		int gridj = _probePositions( i, 1 );
		int index = i + 1; // index of the probe
		int gridi = getProbeIndexY( index );
		int gridj = getProbeIndexY( index );
		fprintf( _fprobes[i], "%5d %.5e %.14e %.14e %.14e %.14e %.14e\n",  
				 state.timestep, state.time,  
				 u(_lev, gridi, gridj), v(_lev, gridi, gridj),
				 state.q(_lev, X, gridi, gridj), state.q(_lev, Y, gridi, gridj),
				 state.omega(_lev, gridi, gridj) ); 
		fflush( _fprobes[i] );
	}
	
	return true;
}

// Update the 2d Array _probePositions to include the gridpoint indices of 
// the currently added probe.  
void OutputProbes::addProbeByIndex( int i, int j ){
	if ( _flagInitialization == true ) {
		cout << "Warning: Addition of probes is allowed only before initialization." << endl;
		exit(1);
	}
	
	if ( ( i < 1  ) || ( j < 1 ) || ( i > _grid.Nx() - 1) || ( j > _grid.Ny() - 1)) {
		cout<< "Warning: invalid probe position:" << 
			" Probes should be located at the inner nodes at the finest grid level." << endl;
		exit(1);
	};
	
	_numProbes += 1;
		
	if ( _numProbes == 1 ) {
		_probePositions.Allocate( 1, _dimen);
		_probePositions( 0, 0 ) = i;
		_probePositions( 0, 1 ) = j; 
	}
	else {
		Array::Array2<int> _probeTemp( _numProbes - 1, _dimen );
		assert( _probeTemp.Size() == _probePositions.Size() );
		_probeTemp = _probePositions;
		
		// increase the size of _probePositions, due to addition of probes
		_probePositions.Deallocate();
		_probePositions.Allocate( _numProbes , _dimen );
		
		for (int dim1 = 0; dim1 < _numProbes - 1; dim1 ++ ) {
			_probePositions( dim1, 0 ) = _probeTemp( dim1, 0 );
			_probePositions( dim1, 1 ) = _probeTemp( dim1, 1 );
		}
		
		// assign the gridpoint indices of the currently added probe
		_probePositions( _numProbes - 1 , 0 ) = i;
		_probePositions( _numProbes - 1 , 1 ) = j;
	}
}

void OutputProbes::addProbeByPosition( double xcord, double ycord ) {
	int i = _grid.getXGridIndex( xcord );
	int j = _grid.getYGridIndex( ycord );
	addProbeByIndex( i, j ); 
}

void OutputProbes::addProbe( int i, int j ) {
	addProbeByIndex( i, j );
}

void OutputProbes::addProbe( double xcord, double ycord ) {
	addProbeByPosition( xcord, ycord ); 
}
	
void OutputProbes::print() {
	cout << "\n-- Probe locations: by grid indices --" << endl;
    for (int dim1 = 0; dim1 < _numProbes ; dim1 ++ ) {
		cout << " " << dim1 + 1  << "-th probe: grid index ( " 
			 << _probePositions(dim1, 0) << ", " << _probePositions(dim1, 1) 
			 << " )" << endl;
	}
	cout << "---------- end of probe list --------" << endl;
}

} // namespace ibpm
