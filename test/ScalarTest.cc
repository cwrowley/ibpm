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

class ScalarTest : public testing::Test {
protected:
    ScalarTest() : 
        _nx(4),
        _ny(8),
        _ngrid(3),
        _grid(_nx, _ny, _ngrid, 2, -1, 3),
        _f(_grid),
        _g(_grid),
        _x(_grid),
        _y(_grid)
    {
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

    // functions to test
    double f(int lev, int i, int j) {
        return lev * 100 + i * 10 + j;
    }

    double g(int lev, int i, int j) {
        return f(lev,i,j) * 1000;
    }
    
    // data
    int _nx;
    int _ny;
    int _ngrid;
    Grid _grid;
    Scalar _f;
    Scalar _g;
    Scalar _x;
    Scalar _y;
};

#define EXPECT_ALL_EQUAL(a,b)                   \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for (int i=1; i<_nx; ++i) {             \
            for (int j=1; j<_ny; ++j) {         \
                EXPECT_DOUBLE_EQ((a), (b));     \
            }                                   \
        }                                       \
    }

TEST_F(ScalarTest, GridSizes ) {
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

TEST_F(ScalarTest, Assignment) {
    int i = _nx/2;
    int j = _ny/2;
    _f(0,i,j) = 73;
    EXPECT_DOUBLE_EQ(_f(0,i,j), 73);
}

TEST_F(ScalarTest, CopyConstructor) {
    Scalar h = _f;
    _f = 0.;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) );
}
    
TEST_F(ScalarTest, CopyFromScalar) {
    Scalar h(_grid);
    h = _f;
    _f = 0.;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) );
}

TEST_F(ScalarTest, CopyFromDouble) {
    double a = 7;
    _f = a;
    EXPECT_ALL_EQUAL( _f(lev,i,j), a );
}

TEST_F(ScalarTest, Elements) {
    EXPECT_ALL_EQUAL( _f(lev,i,j), f(lev,i,j) );
}

TEST_F(ScalarTest, PlusEqualsScalar) {
    _g += _f;
    EXPECT_ALL_EQUAL( _g(lev,i,j), f(lev,i,j) + g(lev,i,j) );
}

TEST_F(ScalarTest, PlusEqualsDouble) {
    _g += 6;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) + 6 );
}

TEST_F(ScalarTest, ScalarPlusScalar) {
    Scalar h = _f + _g;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) + g(lev,i,j) );
}

TEST_F(ScalarTest, ScalarPlusDouble) {
    Scalar h = _f + 4;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) + 4 );
}

TEST_F(ScalarTest, DoublePlusScalar) {
    Scalar h = 4 + _f;
    EXPECT_ALL_EQUAL( h(lev,i,j), f(lev,i,j) + 4 );
}

TEST_F(ScalarTest, MinusEquals) {
    _g -= _f;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) - f(lev,i,j) );
}

TEST_F(ScalarTest, ScalarMinusScalar) {
    Scalar h = _g - _f;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) - f(lev,i,j) );
}

TEST_F(ScalarTest, ScalarMinusDouble) {
    Scalar h = _g - 9;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) - 9 );
}

TEST_F(ScalarTest, DoubleMinusScalar) {
    Scalar h = 9 - _g;
    EXPECT_ALL_EQUAL( h(lev,i,j), 9 - g(lev,i,j) );
}

TEST_F(ScalarTest, TimesEquals) {
    _g *= _f;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) * f(lev,i,j) );
}

TEST_F(ScalarTest, TimesEqualsDouble) {
    _g *= 3;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) * 3 );
}

TEST_F(ScalarTest, ScalarTimesScalar) {
    Scalar h = _g * _f;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) * f(lev,i,j) );
}

TEST_F(ScalarTest, ScalarTimesDouble) {
    Scalar h = _g * 11;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) * 11 );
}

TEST_F(ScalarTest, DoubleTimesScalar) {
    Scalar h = 11 * _g;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) * 11 );
}

TEST_F(ScalarTest, DivEquals) {
    _g /= _f;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) / f(lev,i,j) );
}

TEST_F(ScalarTest, DivEqualsDouble) {
    _g /= 3;
    EXPECT_ALL_EQUAL( _g(lev,i,j), g(lev,i,j) / 3 );
}

TEST_F(ScalarTest, ScalarDivScalar) {
    Scalar h = _g / _f;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) / f(lev,i,j) );
}

TEST_F(ScalarTest, ScalarDivDouble) {
    Scalar h = _g / 11;
    EXPECT_ALL_EQUAL( h(lev,i,j), g(lev,i,j) / 11 );
}

TEST_F(ScalarTest, DoubleDivScalar) {
    Scalar h = 11 / _g;
    EXPECT_ALL_EQUAL( h(lev,i,j), 11 / g(lev,i,j) );
}

TEST_F(ScalarTest, UnaryMinus) {
    Scalar h = -_f;
    EXPECT_ALL_EQUAL( h(lev,i,j), -f(lev,i,j) );
}
    
TEST_F(ScalarTest, SliceAccess) {
    Array2<double> f0 = _f[0];
    Array2<double> f1 = _f[1];
    Array2<double> f2 = _f[2];
    EXPECT_ALL_EQUAL( f0(i,j), _f(0,i,j) );
    EXPECT_ALL_EQUAL( f1(i,j), _f(1,i,j) );
    EXPECT_ALL_EQUAL( f2(i,j), _f(2,i,j) );    
}
    
TEST_F(ScalarTest, SliceAssignment) {
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

    TEST_F(ScalarTest, GetBC) {
        BC bcx(_nx,_ny);
        BC bcy(_nx,_ny);
        // get boundary values on finest grid, from coarsest grid
        int lev=0;
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
    TEST_F(ScalarTest, CoarsifyConstant) {
        double val = 3.;
        _f = val;
        _f.coarsify();
        EXPECT_ALL_EQUAL( val, _f(lev,i,j) );
    }
    
#undef EXPECT_ALL_EQUAL

}
