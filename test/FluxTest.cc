#include "Grid.h"
#include "Flux.h"
#include <gtest/gtest.h>

#define EXPECT_ALL_X_EQUAL(a,b)           \
	for ( int i=0; i<_nx+1; ++i ) {       \
        for ( int j=0; j<_ny; ++j ) {     \
		    EXPECT_DOUBLE_EQ( (a), (b) ); \
        }                                 \
	}

#define EXPECT_ALL_Y_EQUAL(a,b)           \
	for ( int i=0; i<_nx; ++i ) {         \
        for ( int j=0; j<_ny+1; ++j ) {   \
		    EXPECT_DOUBLE_EQ( (a), (b) ); \
        }                                 \
	}

class FluxTest : public testing::Test {
protected:
    FluxTest() :
        _nx(3),
        _ny(5),
        _grid(_nx, _ny, 2, -1, -3),
        _f(_grid),
        _g(_grid)
    {
		for (int i=0; i<_nx+1; ++i) {
			for (int j=0; j<_ny; ++j) {
				_f(X,i,j) = fx(i,j);
				_g(X,i,j) = gx(i,j);
			}
		}
		for (int i=0; i<_nx; ++i) {
			for (int j=0; j<_ny+1; ++j) {
				_f(Y,i,j) = fy(i,j);
				_g(Y,i,j) = gy(i,j);
			}
		}
    }
    
	inline double fx(int i, int j) {
		return i * _nx + j;
	}
	
	inline double fy(int i, int j) {
		return j * _ny - i;
	}

	inline double gx(int i, int j) {
		return i + 3*j;
	}

	inline double gy(int i, int j) {
		return i*j - 10;
	}
	
	inline double f(int i, int j) {
		return 0.5 * i * i * _nx + 2 * j * (i + 1) + cos(j) * (_ny + 1);
	}

	int _nx;
    int _ny;
    Grid _grid;
    Flux _f;
    Flux _g;
};

TEST_F(FluxTest, Assignment) {
	_f(X,1,3) = 73;
	EXPECT_DOUBLE_EQ(_f(X,1,3), 73 );
}

TEST_F(FluxTest, CopyConstructor) {
	Flux h = _f;
	EXPECT_ALL_X_EQUAL( h(X,i,j), fx(i,j) );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), fy(i,j) );
}

TEST_F(FluxTest, CopyFromDouble) {
	double a = 7;
	_f = a;
	EXPECT_ALL_X_EQUAL( _f(X,i,j), a );
	EXPECT_ALL_Y_EQUAL( _f(Y,i,j), a );
}

TEST_F(FluxTest, UnaryMinus) {
	Flux h = -_f;
	EXPECT_ALL_X_EQUAL( h(X,i,j), -fx(i,j) );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), -fy(i,j) );
}

TEST_F(FluxTest, PlusEqualsFlux) {
	_g += _f;
	EXPECT_ALL_X_EQUAL( _g(X,i,j), fx(i,j) + gx(i,j) );
	EXPECT_ALL_Y_EQUAL( _g(Y,i,j), fy(i,j) + gy(i,j) );
}

TEST_F(FluxTest, PlusEqualsDouble) {
	_g += 7;
	EXPECT_ALL_X_EQUAL( _g(X,i,j), gx(i,j) + 7 );
	EXPECT_ALL_Y_EQUAL( _g(Y,i,j), gy(i,j) + 7 );
}

TEST_F(FluxTest, MinusEqualsFlux) {
	_g -= _f;
	EXPECT_ALL_X_EQUAL( _g(X,i,j), gx(i,j) - fx(i,j) );
	EXPECT_ALL_Y_EQUAL( _g(Y,i,j), gy(i,j) - fy(i,j) );
}

TEST_F(FluxTest, MinusEqualsDouble) {
	_g -= 7;
	EXPECT_ALL_X_EQUAL( _g(X,i,j), gx(i,j) - 7 );
	EXPECT_ALL_Y_EQUAL( _g(Y,i,j), gy(i,j) - 7 );
}

TEST_F(FluxTest, FluxPlusFlux) {
	Flux h = _f + _g;
	EXPECT_ALL_X_EQUAL( h(X,i,j), fx(i,j) + gx(i,j) );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), fy(i,j) + gy(i,j) );
}

TEST_F(FluxTest, FluxPlusDouble) {
	Flux h = _f + 7;
	EXPECT_ALL_X_EQUAL( h(X,i,j), fx(i,j) + 7 );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), fy(i,j) + 7 );		
}

TEST_F(FluxTest, DoublePlusFlux) {
	Flux h = 7 + _f;
	EXPECT_ALL_X_EQUAL( h(X,i,j), fx(i,j) + 7 );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), fy(i,j) + 7 );		
}

TEST_F(FluxTest, FluxMinusFlux) {
	Flux h = _f - _g;
	EXPECT_ALL_X_EQUAL( h(X,i,j), fx(i,j) - gx(i,j) );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), fy(i,j) - gy(i,j) );
}

TEST_F(FluxTest, FluxMinusDouble) {
	Flux h = _f - 7;
	EXPECT_ALL_X_EQUAL( h(X,i,j), fx(i,j) - 7 );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), fy(i,j) - 7 );		
}

TEST_F(FluxTest, DoubleMinusFlux) {
	Flux h = 7 - _f;
	EXPECT_ALL_X_EQUAL( h(X,i,j), 7 - fx(i,j) );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), 7 - fy(i,j) );		
}

TEST_F(FluxTest, TimesEqualsDouble) {
	_f *= 3;
	EXPECT_ALL_X_EQUAL( _f(X,i,j), 3 * fx(i,j) );
	EXPECT_ALL_Y_EQUAL( _f(Y,i,j), 3 * fy(i,j) );
}

TEST_F(FluxTest, DivEqualsDouble) {
	_f /= 3.;
	EXPECT_ALL_X_EQUAL( _f(X,i,j), fx(i,j) / 3. );
	EXPECT_ALL_Y_EQUAL( _f(Y,i,j), fy(i,j) / 3. );
}

TEST_F(FluxTest, FluxTimesDouble) {
	double a = 5;
	Flux h = _f * a;
	EXPECT_ALL_X_EQUAL( h(X,i,j), fx(i,j) * a );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), fy(i,j) * a );		
}

TEST_F(FluxTest, DoubleTimesFlux) {
	double a = 5;
	Flux h = a * _f;
	EXPECT_ALL_X_EQUAL( h(X,i,j), fx(i,j) * a );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), fy(i,j) * a );		
}

TEST_F(FluxTest, FluxDivDouble) {
	double a = 5;
	Flux h = _f / a;
	EXPECT_ALL_X_EQUAL( h(X,i,j), fx(i,j) / a );
	EXPECT_ALL_Y_EQUAL( h(Y,i,j), fy(i,j) / a );		
}	

TEST_F(FluxTest, DotProductSymmetric) {
	EXPECT_DOUBLE_EQ( _f.dot(_g), _g.dot(_f) );
}

TEST_F(FluxTest, DotProductDomainArea) {
	Flux h( _grid );
	Flux l( _grid );
	for (int i=0; i<_nx+1; ++i) {
		for (int j=0; j<_ny; ++j) {
			h(X,i,j) = 1;
			l(X,i,j) = 0.2;
		}
	}
	for (int i=0; i<_nx; ++i) {
		for (int j=0; j<_ny+1; ++j) {
			h(Y,i,j) = 0.8;
			l(Y,i,j) = 1;
		}
	}
	double innerarea = (_grid).getArea() * 
						(_nx - 1) *( _ny - 1) / ( _nx * _ny );
	EXPECT_DOUBLE_EQ(h.dot(l), innerarea);
}
	
TEST_F(FluxTest, DotProductValue) {
	double dp = 0;
	for (int i = 1; i < _nx; ++i) {
		for ( int j = 1; j < _ny; ++j) {
			dp += fx(i, j) * gx(i, j);
			dp += fy(i, j) * gy(i, j);
		}			
	}		
	dp *= pow(_grid.getDx(), 2);	
	EXPECT_DOUBLE_EQ( _f.dot(_g), dp ); 
}

// void testIteratorCount() {
//     Flux::iterator i;
//     int n;
// 
//     // X direction
//     // Number of fluxes in X-dir is (nx+1, ny)
//     n = 0;
//     for (i = _f->begin(Flux::X); i != _f->end(Flux::X); ++i) {
//         ++n;
//     }
//     TS_ASSERT_EQUALS(n, (_nx + 1) * _ny );
// 
//     // Y direction
//     // Number of fluxes in Y-dir is (nx, ny+1)
//     n = 0;
//     for (i = _f->begin(Flux::Y); i != _f->end(Flux::Y); ++i) {
//         ++n;
//     }
//     TS_ASSERT_EQUALS(n, _nx * (_ny + 1) );
// }
// 
// void testIteratorAccess() {
//     Flux::iterator iter;
//     int i=0;
//     int j=0;
// 
//     // X dimension
//     for (iter=_f->begin(Flux::X); iter != _f->end(Flux::X); ++iter) {
//         TS_ASSERT_DELTA( *iter, fx(i,j), _delta );
//         ++j;
//         if (j >= _ny) {
//             j -= _ny;
//             ++i;
//         }
//     }
//     
//     // Y dimension
//     // Caution: number of fluxes in Y-dir is ny+1 !
//     i = 0;
//     j = 0;
//     for (iter=_f->begin(Flux::Y); iter != _f->end(Flux::Y); ++iter) {
//         TS_ASSERT_DELTA( *iter, fy(i,j), _delta );
//         ++j;
//         if (j >= (_ny+1)) {
//             j -= (_ny+1);
//             ++i;
//         }
//     }
// }
// 
// void testIteratorAssignment() {
//     Flux::iterator iter;
//     int i=0;
//     int j=0;
// 
//     // X dimension
//     for (iter=_f->begin(Flux::X); iter != _f->end(Flux::X); ++iter) {
//         *iter = 2 * fx(i,j) + 3;
//         ++j;
//         if (j >= _ny) {
//             j -= _ny;
//             ++i;
//         }
//     }
//     ASSERT_ALL_X_EQUAL( _f->x(i,j), 2 * fx(i,j) + 3 );
//     ASSERT_ALL_Y_EQUAL( _f->y(i,j), fy(i,j) );
// 
//     // Y dimension
//     // Caution: number of fluxes in Y-dir is ny+1 !
//     i = 0;
//     j = 0;
//     for (iter=_g->begin(Flux::Y); iter != _g->end(Flux::Y); ++iter) {
//         *iter = 3 * gy(i,j) + 2;
//         ++j;
//         if (j >= (_ny+1)) {
//             j -= (_ny+1);
//             ++i;
//         }
//     }
//     ASSERT_ALL_X_EQUAL( _g->x(i,j), gx(i,j) );
//     ASSERT_ALL_Y_EQUAL( _g->y(i,j), 3 * gy(i,j) + 2 );
// 
// }


#undef EXPECT_ALL_X_EQUAL
#undef EXPECT_ALL_Y_EQUAL
