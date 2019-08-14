#ifndef _STATEVECTOR_H_
#define _STATEVECTOR_H_

#include "State.h"
#include "Grid.h"

namespace ibpm {

/*!
    \file StateVector.h
    \class StateVector

    \brief Structure for grouping state variables as a vector

    \author Daniel Floryan
    \author $LastChangedBy$
    \date  16 Jan 2019
    \date $LastChangedDate$
    \version $Revision$
*/

class StateVector {
public:
    /// Default constructor: do not allocate memory
    StateVector();

    StateVector( const StateVector& v);

    StateVector( const State& y);
    	  		
    StateVector( const Grid& grid, int numPoints );

    /// \brief Instantiate a state by reading data from the specified file
    StateVector( string filename );

    ~StateVector();
    
    /// \brief Allocate memory, with the specified Grid and number of
    /// boundary points
    void resize( const Grid& grid, int numPoints );

    /// \brief Save the State to a file (e.g. as a restart file)
    /// Return true if successful
    /*  WARNING:  At this point, the xshift and yshift parameters are not saved
     and are not checked for compatibility when loading.  Caution should be
     taken when working with shifted grids.  This approach was taken to pre-
     serve backwards compatibility with previously saved binary files.  In 
     the future perhaps using HDF5 would prevent such problems.
     */
    bool save(std::string filename) const;

    /// \brief Load the state from a file (e.g. as a restart file)
    /// Return true if successful
    bool load(const std::string& filename);

    /// Copy assignment
    inline StateVector& operator=(const StateVector& v) {
        x.q = v.x.q;
        x.omega = v.x.omega;
        x.f = v.x.f;
        return *this;
    }

    /// Copy assignment from double
    inline StateVector& operator=(double a) {
        x.q = a;
        x.omega = a;
        x.f = a;
        return *this;
    }

    /// u += v
    inline StateVector& operator+=(const StateVector& v) {
        x.q += v.x.q;
        x.omega += v.x.omega;
        x.f += v.x.f;
        return *this;
    }

    /// u += a
    inline StateVector& operator+=(double a) {
        x.q += a;
        x.omega += a;
        x.f += a;
        return *this;
    }

    /// u -= v
    inline StateVector& operator-=(const StateVector& v) {
        x.q -= v.x.q;
        x.omega -= v.x.omega;
        x.f -= v.x.f;
        return *this;
    }

    /// u -= a
    inline StateVector& operator-=(double a) {
        x.q -= a;
        x.omega -= a;
        x.f -= a;
        return *this;
    }

    /// u *= a
    inline StateVector& operator*=(double a) {
        x.q *= a;
        x.omega *= a;
        x.f *= a;
        return *this;
    }

    /// u /= a
    inline StateVector& operator/=(double a) {
        x.q /= a;
        x.omega /= a;
        x.f /= a;
        return *this;
    }

    /// u + v
    inline StateVector operator+(const StateVector& v) {
        StateVector w = *this;
        w += v;
        return w;
    }

    /// u + a
    inline StateVector operator+(double a) {
        StateVector w = *this;
        w += a;
        return w;
    }
    
    /// u - v
    inline StateVector operator-(const StateVector& v) {
        StateVector w = *this;
        w -= v;
        return w;
    }

    /// u - a
    inline StateVector operator-(double a) {
        StateVector w = *this;
        w -= a;
        return w;
    }
    
    /// u * a
    inline StateVector operator*(double a) {
        StateVector w = *this;
        w *= a;
        return w;
    }
    
    /// u / a
    inline StateVector operator/(double a) {
        StateVector w = *this;
        w /= a;
        return w;
    }
			    
    // public data
    State x;
};  // class StateVector

/// -v
inline StateVector operator-(const StateVector& v) {
    StateVector w = v;
    w *= -1;
    return w;
}

/// u + v
inline StateVector operator+(const StateVector& u, const StateVector& v) {
    StateVector w = u;
    w += v;
    return w;
};

/// a + v
inline StateVector operator+(double a, const StateVector& v) {
    StateVector w = v;
    w += a;
    return w;
};

/// u - v
inline StateVector operator-(const StateVector& u, const StateVector& v) {
    StateVector w = u;
    w -= v;
    return w;
};
    
/// a - v
inline StateVector operator-(double a, const StateVector& v) {
    StateVector w = -v;
    w += a;
    return w;
};

/// a * v
inline StateVector operator*(double a, const StateVector& v) {
    StateVector w = v;
    w *= a;
    return w;
}

/// v * a
inline StateVector operator*(const StateVector& v, double a) {
    StateVector w = v;
    w *= a;
    return w;
}

/// v / a
inline StateVector operator/(const StateVector& v, double a) {
    StateVector w = v;
    w /= a;
    return w;
}

} // namespace ibpm

#endif /* _STATEVECTOR_H_ */
