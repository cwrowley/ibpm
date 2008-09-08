// bininfo - utility for printing header information from a restart file
//
// Clancy Rowley
// 7 Sep 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include <iostream>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
    }
    
    State x;

    // Read in a restart file
    if ( ! x.load( argv[1] ) ) {
        cerr << "Error reading file " << argv[1] << endl;
    }
    else {
        cout << "Timestep " << x.timestep << endl;
        cout << "time = " << x.time << endl;
        // Write out the Grid information
        cout << "nx = " << x.gamma.getNx() << endl
             << "ny = " << x.gamma.getNy() << endl
             << "dx = " << x.gamma.getDx() << endl
             << "x0 = " << x.gamma.getXEdge(0) << endl
             << "y0 = " << x.gamma.getYEdge(0) << endl
             << "numPoints = " << x.f.getNumPoints() << endl;
    }
    
    return 0;
}