#include "Grid.h"
#include "Flux.h"
#include <math.h>
#include <gtest/gtest.h>

using namespace ibpm;

#define EXPECT_ALL_X_EQUAL(a,b)           \
    for ( int i=0; i<_nx+1; ++i ) {       \
        for ( int j=0; j<_ny; ++j ) {     \
            EXPECT_DOUBLE_EQ( (a), (b) ); \
        }                                 \
    }

#define EXPECT_ALL_X_NEAR(a,b,tol)        \
    for ( int i=0; i<_nx+1; ++i ) {       \
        for ( int j=0; j<_ny; ++j ) {     \
            EXPECT_NEAR( (a), (b), tol ); \
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
    _f = 0.;
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

TEST_F(FluxTest, FluxCoordinates) {
    EXPECT_ALL_X_EQUAL( _f.x(X,i), _grid.getXEdge(i) );
    EXPECT_ALL_X_EQUAL( _f.y(X,j), _grid.getYCenter(j) );
    EXPECT_ALL_Y_EQUAL( _f.x(Y,i), _grid.getXCenter(i) );
    EXPECT_ALL_Y_EQUAL( _f.y(Y,j), _grid.getYEdge(j) );
}

TEST_F(FluxTest, IndexCount) {
    Flux::index ind;
    int count=0;
    for (ind = _f.begin(); ind != _f.end(); ++ind) {
        ++count;
    }
    EXPECT_EQ(count, 2*_nx*_ny + _nx + _ny);

    count=0;
    for (ind = _f.begin(X); ind != _f.end(X); ++ind) {
        ++count;
    }
    EXPECT_EQ(count, _nx*_ny + _ny);
    
    count=0;
    for (ind = _f.begin(Y); ind != _f.end(Y); ++ind) {
        ++count;
    }
    EXPECT_EQ(count, _nx*_ny + _nx);
}

TEST_F(FluxTest, IndexReference) {
    Flux::index ind;

    // loop over all values
    for (ind = _f.begin(); ind != _f.end(); ++ind) {
        _f(ind) = 35;
    }
    EXPECT_ALL_X_EQUAL( _f(X,i,j), 35 );
    EXPECT_ALL_Y_EQUAL( _f(Y,i,j), 35 );
    for (ind = _f.begin(); ind != _f.end(); ++ind) {
        EXPECT_DOUBLE_EQ( _f(ind), 35 );
    }
    
    // loop over X alone
    for (ind = _f.begin(X); ind != _f.end(X); ++ind) {
        _f(ind) = 53;
    }
    EXPECT_ALL_X_EQUAL( _f(X,i,j), 53 );
    EXPECT_ALL_Y_EQUAL( _f(Y,i,j), 35 );

    // loop over Y alone
    for (ind = _f.begin(Y); ind != _f.end(Y); ++ind) {
        _f(ind) = 350;
    }
    EXPECT_ALL_X_EQUAL( _f(X,i,j), 53 );
    EXPECT_ALL_Y_EQUAL( _f(Y,i,j), 350 );

}

TEST_F(FluxTest, GetIndex) {
    EXPECT_ALL_X_EQUAL( _f(_f.getIndex(X,i,j)), _f(X,i,j) );
    EXPECT_ALL_Y_EQUAL( _f(_f.getIndex(Y,i,j)), _f(Y,i,j) );
}

TEST_F(FluxTest, GridValues) {
    EXPECT_ALL_X_EQUAL( _f.x(_f.getIndex(X,i,j)), _f.x(X,i) );
    EXPECT_ALL_Y_EQUAL( _f.x(_f.getIndex(Y,i,j)), _f.x(Y,i) );    
    EXPECT_ALL_X_EQUAL( _f.y(_f.getIndex(X,i,j)), _f.y(X,j) );
    EXPECT_ALL_Y_EQUAL( _f.y(_f.getIndex(Y,i,j)), _f.y(Y,j) );    
}

TEST_F(FluxTest, UniformFlow) {
    double mag = 3.;
    double angle = 0.;
    Flux q = Flux::UniformFlow( _grid, mag, angle );
    EXPECT_ALL_X_EQUAL( q(X,i,j), mag * _grid.getDx() );
    EXPECT_ALL_Y_EQUAL( q(Y,i,j), 0. );

    double pi = 4. * atan(1.);
    angle = pi/2.;
    q = Flux::UniformFlow( _grid, mag, angle );
    EXPECT_ALL_X_NEAR( q(X,i,j), 0., 1e-15 );
    EXPECT_ALL_Y_EQUAL( q(Y,i,j), mag * _grid.getDx() );
}

#undef EXPECT_ALL_X_EQUAL
#undef EXPECT_ALL_Y_EQUAL
