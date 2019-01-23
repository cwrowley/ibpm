// StateVector.cc
//
// Description:
// Implementation of the StateVector class
//
// Author(s):
// Daniel Floryan
//
// Date: 16 Jan 2019
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "StateVector.h"
#include <string>
#include <stdio.h>

using namespace ibpm;

namespace ibpm {

StateVector::StateVector() {}

StateVector::StateVector(const StateVector& v) {
    resize( v.x.omega.getGrid(), v.x.f.getNumPoints() );
    x.q = v.x.q;
    x.omega = v.x.omega;
    x.f = v.x.f;
}

StateVector::StateVector(const State& y) {
    resize( y.omega.getGrid(), y.f.getNumPoints() );
    x.q = y.q;
    x.omega = y.omega;
    x.f = y.f;
}

StateVector::StateVector(const Grid& grid, int numPoints ) {
    resize( grid, numPoints );
}

void StateVector::resize( const Grid& grid, int numPoints ) {
    x.q.resize( grid );
    x.omega.resize( grid );
    x.f.resize( numPoints );
}

StateVector::StateVector( string filename ) {
    load( filename );
}

StateVector::~StateVector() {}

bool StateVector::load(const std::string& filename) {
    return x.load(filename);
}

bool StateVector::save(std::string filename) const {
    return x.save(filename);
}

} // namespace ibpm
