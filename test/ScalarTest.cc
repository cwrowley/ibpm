#include "Grid.h"
#include "Scalar.h"
#include <gtest/gtest.h>

namespace {

class ScalarTest : public testing::Test {
protected:
    ScalarTest() : 
        _nx(3),
        _ny(5),
        _grid(_nx, _ny, 2, -1, 3),
        _f(_grid),
        _g(_grid)
    {
  		for (int i=0; i<_nx+1; ++i) {
  			for (int j=0; j<_ny+1; ++j) {
  				_f(i,j) = f(i,j);
  				_g(i,j) = g(i,j);
  			}
  		}
  		
    }

	// functions to test
	double f(int i, int j) {
		return 0.5 * i * i * _nx + 2 * j * (i + 1) + cos(j) * (_ny + 1);
	}

	double g(int i, int j) {
		return 3*i + j;
	}
	
	double fx(int i, int j) {
		return 3 * i * i * _nx + 5 * (j +1 )* (1 + i + 0.7 * j) * cos(j);
	}
	
	double fy(int i, int j) {
		return  i * (i + j ) *3  + (_ny + cos(i)) * (j + 1 ) *2;
	}

    // data
	int _nx;
	int _ny;
	Grid _grid;
	Scalar _f;
	Scalar _g;
};

#define EXPECT_ALL_EQUAL(a,b)               \
	for (int i=0; i<_nx+1; ++i) {           \
		for (int j=0; j<_ny+1; ++j) {       \
			EXPECT_DOUBLE_EQ((a), (b));     \
		}                                   \
	}

TEST_F(ScalarTest, Assignment) {
	int i = _nx/2;
	int j = _ny/2;
	_f(i,j) = 73;
	EXPECT_DOUBLE_EQ(_f(i,j), 73);
}

TEST_F(ScalarTest, CopyConstructor) {
	Scalar h = _f;
    _f = 0.;
	EXPECT_ALL_EQUAL( h(i,j), f(i,j) );
}

TEST_F(ScalarTest, CopyFromDouble) {
	double a = 7;
	_f = a;
	EXPECT_ALL_EQUAL( _f(i,j), a );
}

TEST_F(ScalarTest, Elements) {
	EXPECT_ALL_EQUAL( _f(i,j), f(i,j) );
}

TEST_F(ScalarTest, PlusEqualsScalar) {
	_g += _f;
	EXPECT_ALL_EQUAL( _g(i,j), f(i,j) + g(i,j) );
}

TEST_F(ScalarTest, PlusEqualsDouble) {
	_g += 6;
	EXPECT_ALL_EQUAL( _g(i,j), g(i,j) + 6 );
}

TEST_F(ScalarTest, ScalarPlusScalar) {
	Scalar h = _f + _g;
	EXPECT_ALL_EQUAL( h(i,j), f(i,j) + g(i,j) );
}

TEST_F(ScalarTest, ScalarPlusDouble) {
	Scalar h = _f + 4;
	EXPECT_ALL_EQUAL( h(i,j), f(i,j) + 4 );
}

TEST_F(ScalarTest, DoublePlusScalar) {
	Scalar h = 4 + _f;
	EXPECT_ALL_EQUAL( h(i,j), f(i,j) + 4 );
}

TEST_F(ScalarTest, MinusEquals) {
	_g -= _f;
	EXPECT_ALL_EQUAL( _g(i,j), g(i,j) - f(i,j) );
}

TEST_F(ScalarTest, ScalarMinusScalar) {
	Scalar h = _g - _f;
	EXPECT_ALL_EQUAL( h(i,j), g(i,j) - f(i,j) );
}

TEST_F(ScalarTest, ScalarMinusDouble) {
	Scalar h = _g - 9;
	EXPECT_ALL_EQUAL( h(i,j), g(i,j) - 9 );
}

TEST_F(ScalarTest, DoubleMinusScalar) {
	Scalar h = 9 - _g;
	EXPECT_ALL_EQUAL( h(i,j), 9 - g(i,j) );
}

TEST_F(ScalarTest, TimesEquals) {
	_g *= _f;
	EXPECT_ALL_EQUAL( _g(i,j), g(i,j) * f(i,j) );
}

TEST_F(ScalarTest, TimesEqualsDouble) {
	_g *= 3;
	EXPECT_ALL_EQUAL( _g(i,j), g(i,j) * 3 );
}

TEST_F(ScalarTest, ScalarTimesScalar) {
	Scalar h = _g * _f;
	EXPECT_ALL_EQUAL( h(i,j), g(i,j) * f(i,j) );
}

TEST_F(ScalarTest, ScalarTimesDouble) {
	Scalar h = _g * 11;
	EXPECT_ALL_EQUAL( h(i,j), g(i,j) * 11 );
}

TEST_F(ScalarTest, DoubleTimesScalar) {
	Scalar h = 11 * _g;
	EXPECT_ALL_EQUAL( h(i,j), g(i,j) * 11 );
}

TEST_F(ScalarTest, DivEquals) {
	_g /= _f;
	EXPECT_ALL_EQUAL( _g(i,j), g(i,j) / f(i,j) );
}

TEST_F(ScalarTest, DivEqualsDouble) {
	_g /= 3;
	EXPECT_ALL_EQUAL( _g(i,j), g(i,j) / 3 );
}

TEST_F(ScalarTest, ScalarDivScalar) {
	Scalar h = _g / _f;
	EXPECT_ALL_EQUAL( h(i,j), g(i,j) / f(i,j) );
}

TEST_F(ScalarTest, ScalarDivDouble) {
	Scalar h = _g / 11;
	EXPECT_ALL_EQUAL( h(i,j), g(i,j) / 11 );
}

TEST_F(ScalarTest, DoubleDivScalar) {
	Scalar h = 11 / _g;
	EXPECT_ALL_EQUAL( h(i,j), 11 / g(i,j) );
}

TEST_F(ScalarTest, UnaryMinus) {
	Scalar h = -_f;
	EXPECT_ALL_EQUAL( h(i,j), -f(i,j) );
}
	

//    void testIteratorCount() {
//        Scalar::iterator i;
//        int n=0;
//        for (i = (*_f).begin(); i != (*_f).end(); ++i) {
//            ++n;
//        }
//        TS_ASSERT_EQUALS(n, _nx * _ny);
//    }
//
//    void testIteratorAccess() {
//        Scalar::iterator iter;
//        int i=0;
//        int j=0;
//        for (iter = (*_f).begin(); iter != (*_f).end(); ++iter) {
//            TS_ASSERT_DELTA( *iter, f(i,j), _delta );
//            ++j;
//            if (j >= _ny) {
//                j -= _ny;
//                ++i;
//            }
//        }
//    }
//
//    void testIteratorAssignment() {
//        Scalar::iterator iter;
//        int i=0;
//        int j=0;
//        for (iter = (*_f).begin(); iter != (*_f).end(); ++iter) {
//            *iter = f(i,j) * 2 + 3;
//            ++j;
//            if (j >= _ny) {
//                j -= _ny;
//                ++i;
//            }
//        }
//        EXPECT_ALL_EQUAL( (*_f)(i,j), f(i,j) * 2 + 3 );
//    }

#undef EXPECT_ALL_EQUAL

}
