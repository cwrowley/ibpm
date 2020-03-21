// header file for storing integration schemes

#include "Array.h"
#include <string>

using std::string;
using std::cout;

class Scheme {
public:
	enum SchemeType { EULER, AB2, RK3, RK3b };
	
    Scheme( SchemeType scheme ) { 
        if ( scheme == EULER ) {
            _coeff.Allocate(1,3);
            _coeff(0,0) = 1.;		// an
            _coeff(0,1) = 0.;		// bn
			_coeff(0,2) = 1.;		// cn
			_name = "Explicit Euler";
        }
        else if ( scheme == AB2 ) {
            _coeff.Allocate(1,3);
            _coeff(0,0) = 3./2.;	// an
            _coeff(0,1) = -1./2.;	// bn
			_coeff(0,2) = 1.;		// cn
			_name = "Adams Bashforth";
        }
        else if ( scheme == RK3 ) {
            _coeff.Allocate(3,3);
            _coeff(0,0) = 8./15.;	// an
            _coeff(1,0) = 5./12.;
            _coeff(2,0) = 3./4.;
            _coeff(0,1) = 0.;		// bn
            _coeff(1,1) = -17./60.;
            _coeff(2,1) = -5./12.;
			_coeff(0,2) = 8./15.;	// cn
			_coeff(1,2) = 2./3.;
			_coeff(2,2) = 1.;
			_name = "3rd-order Runge Kutta (3-step)";
        }
        else if ( scheme == RK3b ) {
            _coeff.Allocate(4,3);
            _coeff(0,0) = 8./17.;	// an
            _coeff(1,0) = 17./60.;
            _coeff(2,0) = 5./12.;
            _coeff(3,0) = 3./4.;
            _coeff(0,1) = 0.;		// bn
            _coeff(1,1) = -15./68.;
            _coeff(2,1) = -17./60.;
            _coeff(3,1) = -5./12.;
			_coeff(0,2) = 8./17.;	// cn
			_coeff(1,2) = 8./15.;
			_coeff(2,2) = 2./3.;
            _coeff(3,2) = 1.;
			_name = "3rd-order Runge Kutta (4-step)";
        }
		else {
			cout << endl << "ERROR: unrecognized solver: " << scheme << endl << endl;
			exit(1);
		}
    }
	
    inline double an( int i ) { return _coeff(i,0); }
    inline double bn( int i ) { return _coeff(i,1); }
	inline double cn( int i ) { return _coeff(i,2); }
    inline int nsteps() { return _coeff.Nx(); }
	inline string name() { return _name; }
	
private:
	Array::Array2<double> _coeff;
	string _name;
	
};


