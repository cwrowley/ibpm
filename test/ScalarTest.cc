#include "Array.h"
#include "BC.h"
#include "Grid.h"
#include "Scalar.h"
#include <gtest/gtest.h>
#include <iostream>

using namespace std;
using namespace Array;
using namespace ibpm;

namespace {
	
// Test an x, y shift for complete shifting, shifting by one grid pt, and no shift
// The shift values must be global for the parameterized tests
// According _nx and _ny must also be global	
const int _nx(8);
const int _ny(12);
	
//const double _xShiftVal[] = {-1., (double) -4/_nx, 0., (double) 4/_nx, 1.};	
//const double _yShiftVal[] = {-1., (double) -4/_ny, (double) 4/_ny, 1.};	

const double _xShiftVal[] = {(double) -4/_nx, 0., (double) 4/_nx};	
const double _yShiftVal[] = {(double) -4/_ny, (double) 4/_ny};	

//const double _xShiftVal[] = {(double) -1., 1.};

class ScalarTestX : public testing::TestWithParam<double> {
protected:
    ScalarTestX() : 
        _ngrid(3),
        _grid(_nx, _ny, _ngrid, 2, -1, 3 ),
        _f(_grid),
		_g(_grid),
		_x(_grid),
		_y(_grid)
    { 
        setScalars();
    }
    
	virtual ~ScalarTestX() {}
	
    // functions to test
    double f(int lev, int i, int j) {
        return lev * 100 + i * 10 + j;
    }

    double g(int lev, int i, int j) {
        return f(lev,i,j) * 1000;
    }
	
	// functions
	void TestBC();
	void TestCoarsify();
	void setScalars();
    void resizeScalars();
	
    // data
    int _ngrid;
    Grid _grid;
    Scalar _f;
    Scalar _g;
    Scalar _x;
    Scalar _y;
};
	
	
class ScalarTestY : public ScalarTestX {
protected:
	ScalarTestY() { }
};	
	
#define EXPECT_ALL_EQUAL(a,b)                   \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for (int i=1; i<_nx; ++i) {             \
            for (int j=1; j<_ny; ++j) {         \
                EXPECT_DOUBLE_EQ((a), (b));     \
            }                                   \
        }                                       \
    }

TEST_F(ScalarTestX, GridSizes ) {
    double length = 5.;
    double xOffset = 0;
    double yOffset = 0;
    
    Grid grid( _nx, _ny, _ngrid, length, xOffset, yOffset );
    Scalar u( grid );
    
    EXPECT_EQ( _nx, u.Nx() );
    EXPECT_EQ( _ny, u.Ny() );
    EXPECT_EQ( _ngrid, u.Ngrid() );
    EXPECT_DOUBLE_EQ( grid.Dx(), u.Dx() );
    EXPECT_ALL_EQUAL( grid.getXEdge(lev,i), u.getXEdge(lev,i) );
    EXPECT_ALL_EQUAL( grid.getYEdge(lev,j), u.getYEdge(lev,j) );

}

TEST_F(ScalarTestX, Assignment) {
    int i = _nx/2;
    int j = _ny/2;
    _f(0,i,j) = 73;
    EXPECT_DOUBLE_EQ(_f(0,i,j), 73);
}

TEST_F(ScalarTestX, CopyConstructor) {
    Scalar h = _f;
    _f = 0.;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) );
}
    
TEST_F(ScalarTestX, CopyFromScalar) {
    Scalar h(_grid);
    h = _f;
    _f = 0.;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) );
}

TEST_F(ScalarTestX, CopyFromDouble) {
    double a = 7;
    _f = a;
    EXPECT_ALL_EQUAL( _f(lev,i,j), a );
}

TEST_F(ScalarTestX, Elements) {
    EXPECT_ALL_EQUAL( _f(lev,i,j), f(lev,i,j) );
}

TEST_F(ScalarTestX, PlusEqualsScalar) {
    _g += _f;
    EXPECT_ALL_EQUAL( _g(lev,i,j), f(lev,i,j) + g(lev,i,j) );
}

TEST_F(ScalarTestX, PlusEqualsDouble) {
    _g += 6;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) + 6 );
}

TEST_F(ScalarTestX, ScalarPlusScalar) {
    Scalar h = _f + _g;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) + g(lev,i,j) );
}

TEST_F(ScalarTestX, ScalarPlusDouble) {
    Scalar h = _f + 4;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) + 4 );
}

TEST_F(ScalarTestX, DoublePlusScalar) {
    Scalar h = 4 + _f;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) + 4 );
}

TEST_F(ScalarTestX, MinusEquals) {
    _g -= _f;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) - f(lev,i,j) );
}

TEST_F(ScalarTestX, ScalarMinusScalar) {
    Scalar h = _g - _f;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) - f(lev,i,j) );
}

TEST_F(ScalarTestX, ScalarMinusDouble) {
    Scalar h = _g - 9;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) - 9 );
}

TEST_F(ScalarTestX, DoubleMinusScalar) {
    Scalar h = 9 - _g;
    EXPECT_ALL_EQUAL( h(lev,i,j), 9 - g(lev,i,j) );
}

TEST_F(ScalarTestX, TimesEquals) {
    _g *= _f;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) * f(lev,i,j) );
}

TEST_F(ScalarTestX, TimesEqualsDouble) {
    _g *= 3;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) * 3 );
}

TEST_F(ScalarTestX, ScalarTimesScalar) {
    Scalar h = _g * _f;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) * f(lev,i,j) );
}

TEST_F(ScalarTestX, ScalarTimesDouble) {
    Scalar h = _g * 11;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) * 11 );
}

TEST_F(ScalarTestX, DoubleTimesScalar) {
    Scalar h = 11 * _g;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) * 11 );
}

TEST_F(ScalarTestX, DivEquals) {
    _g /= _f;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) / f(lev,i,j) );
}

TEST_F(ScalarTestX, DivEqualsDouble) {
    _g /= 3;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) / 3 );
}

TEST_F(ScalarTestX, ScalarDivScalar) {
    Scalar h = _g / _f;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) / f(lev,i,j) );
}

TEST_F(ScalarTestX, ScalarDivDouble) {
    Scalar h = _g / 11;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) / 11 );
}

TEST_F(ScalarTestX, DoubleDivScalar) {
    Scalar h = 11 / _g;
    EXPECT_ALL_EQUAL( h(lev,i,j), 11 / g(lev,i,j) );
}

TEST_F(ScalarTestX, UnaryMinus) {
    Scalar h = -_f;
    EXPECT_ALL_EQUAL( h(lev,i,j), -f(lev,i,j) );
}
    
TEST_F(ScalarTestX, SliceAccess) {
    Array2<double> f0 = _f[0];
    Array2<double> f1 = _f[1];
    Array2<double> f2 = _f[2];
    EXPECT_ALL_EQUAL( f0(i,j), _f(0,i,j) );
    EXPECT_ALL_EQUAL( f1(i,j), _f(1,i,j) );
    EXPECT_ALL_EQUAL( f2(i,j), _f(2,i,j) );    
}
    
TEST_F(ScalarTestX, SliceAssignment) {
    Array2<double> f0 = _f[0];
    Array2<double> f1 = _f[1];
    Array2<double> f2 = _f[2];
    for (int i=1; i<_nx; ++i) {
        for (int j=1; j<_ny; ++j) {
            f0(i,j) = g(0,i,j);
            f1(i,j) = g(1,i,j);
            f2(i,j) = g(2,i,j);
        }
    }
    // check to see the original array _f was changed
    EXPECT_ALL_EQUAL( g(lev,i,j), _f(lev,i,j) );
}

TEST_P(ScalarTestX, GetBC) {
	_grid.setXShift( GetParam() );
	resizeScalars();
	setScalars();
	TestBC();
}

TEST_P(ScalarTestX, CoarsifyConstant) {
	_grid.setXShift( GetParam() );
	resizeScalars();
	setScalars();
	TestCoarsify();
}

INSTANTIATE_TEST_CASE_P(
	xShiftTests, ScalarTestX, ::testing::ValuesIn(_xShiftVal) 
);	

TEST_P(ScalarTestY, GetBC) {
	_grid.setYShift( GetParam() );
	resizeScalars();
	setScalars();
    TestBC();
}
	
TEST_P(ScalarTestY, CoarsifyConstant) {
	TestCoarsify();
}

INSTANTIATE_TEST_CASE_P(
	yShiftTests, ScalarTestY, ::testing::ValuesIn(_yShiftVal) 
);		

void ScalarTestX::TestBC() {
	BC bcx(_nx,_ny);
	BC bcy(_nx,_ny);
	// get boundary values on finest grid, from coarsest grid
	int lev=1;
	_x.getBC( lev, bcx );
	_y.getBC( lev, bcy );

	#ifdef DEBUG
	_x.print();
	cout << "Boundary of finest grid:" << endl;

    cout << "top: ";
	for (int i=0; i<=_nx; ++i) {
		cout << bcx.top(i) << " ";
	}
	cout << endl << "bottom: ";
	for (int i=0; i<=_nx; ++i) {
		cout << bcx.bottom(i) << " ";
	}
	cout << endl << "left: ";
	for (int j=0; j<=_ny; ++j) {
		cout << bcx.left(j) << " ";
	}
	cout << endl << "right: ";
	for (int j=0; j<=_ny; ++j) {
		cout << bcx.right(j) << " ";
	}
	#endif // DEBUG
	
	for (int i=0; i<=_nx; ++i) {
		EXPECT_DOUBLE_EQ( _grid.getXEdge(lev,i), bcx.top(i) );
		EXPECT_DOUBLE_EQ( _grid.getXEdge(lev,i), bcx.bottom(i) );
		EXPECT_DOUBLE_EQ( _grid.getYEdge(lev,_ny), bcy.top(i) );
		EXPECT_DOUBLE_EQ( _grid.getYEdge(lev,0), bcy.bottom(i) );
	}
	for (int j=0; j<=_ny; ++j) {
		EXPECT_DOUBLE_EQ( _grid.getXEdge(lev,0), bcx.left(j) );
		EXPECT_DOUBLE_EQ( _grid.getXEdge(lev,_nx), bcx.right(j) );
		EXPECT_DOUBLE_EQ( _grid.getYEdge(lev,j), bcy.left(j) );
		EXPECT_DOUBLE_EQ( _grid.getYEdge(lev,j), bcy.right(j) );
	}
}
	
// If f is constant in space, then coarsify(f) should equal f
void ScalarTestX::TestCoarsify() {
	double val = 3.;
	_f = val;
	_f.coarsify();
	EXPECT_ALL_EQUAL( val, _f(lev,i,j) );
}

// Resize scalar fields if the grid properties are changed
void ScalarTestX::resizeScalars( ) {
	_f.resize(_grid);
	_g.resize(_grid);
	_x.resize(_grid);
	_y.resize(_grid);
}
	
// Recalculate values for test scalars
void ScalarTestX::setScalars( ) {
	for (int lev=0; lev < _ngrid; ++lev) {
		for (int i=1; i<_nx; ++i) {
			for (int j=1; j<_ny; ++j) {
				_f(lev,i,j) = f(lev,i,j);
				_g(lev,i,j) = g(lev,i,j);
				_x(lev,i,j) = _grid.getXEdge(lev,i);
				_y(lev,i,j) = _grid.getYEdge(lev,j);
			}
		}
	}    
}
	
#undef EXPECT_ALL_EQUAL

}
