#include "Grid.h"
#include "Flux.h"
#include <math.h>
#include <gtest/gtest.h>

using namespace ibpm;

#define EXPECT_ALL_X_EQUAL(a,b)                 \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for ( int i=0; i<_nx+1; ++i ) {         \
            for ( int j=0; j<_ny; ++j ) {       \
                EXPECT_DOUBLE_EQ( (a), (b) );   \
            }                                   \
        }                                       \
    }

#define EXPECT_ALL_X_NEAR(a,b,tol)                  \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for ( int i=0; i<_nx+1; ++i ) {         \
            for ( int j=0; j<_ny; ++j ) {       \
                EXPECT_NEAR( (a), (b), tol );   \
            }                                   \
        }                                       \
}

#define EXPECT_ALL_Y_EQUAL(a,b)                 \
    for (int lev=0; lev<_ngrid; ++lev) {        \
        for ( int i=0; i<_nx; ++i ) {           \
            for ( int j=0; j<_ny+1; ++j ) {     \
                EXPECT_DOUBLE_EQ( (a), (b) );   \
            }                                   \
        }                                       \
    }

class FluxTest : public testing::Test {
protected:
    FluxTest() :
        _nx(4),
        _ny(8),
        _ngrid(3),
        _grid(_nx, _ny, _ngrid, 2, -1, -3),
        _f(_grid),
        _g(_grid)
    {
        for (int lev=0; lev<_ngrid; ++lev) {
            for (int i=0; i<_nx+1; ++i) {
                for (int j=0; j<_ny; ++j) {
                    _f(lev,X,i,j) = fx(lev,i,j);
                    _g(lev,X,i,j) = gx(lev,i,j);
                }
            }
        }
        for (int lev=0; lev<_ngrid; ++lev) {
            for (int i=0; i<_nx; ++i) {
                for (int j=0; j<_ny+1; ++j) {
                    _f(lev,Y,i,j) = fy(lev,i,j);
                    _g(lev,Y,i,j) = gy(lev,i,j);
                }
            }
        }
    }
    
    inline double fx(int lev, int i, int j) {
        return (i * 10 + j) * (lev+1);
    }
    
    inline double fy(int lev, int i, int j) {
        return (i * 1000 - 100 * j) * (lev+1);
    }

    inline double gx(int lev, int i, int j) {
        return (i * 0.01 + j * 0.1) * (lev+1);
    }

    inline double gy(int lev, int i, int j) {
        return (i * 0.0001 + j * 0.001 ) * (lev+1);
    }
    
    int _nx;
    int _ny;
    int _ngrid;
    Grid _grid;
    Flux _f;
    Flux _g;
};

TEST_F(FluxTest, Assignment) {
    _f(0,X,1,3) = 73;
    EXPECT_DOUBLE_EQ(_f(0,X,1,3), 73 );
}

TEST_F(FluxTest, GridSizes ) {
    int nx = _nx;
    int ny = _ny;
    int ngrid = _ngrid;
    double length = 5.;
    double xOffset = 0;
    double yOffset = 0;
    
    Grid grid( nx, ny, ngrid, length, xOffset, yOffset );
    Flux q( grid );
    
    EXPECT_EQ( nx, q.Nx() );
    EXPECT_EQ( ny, q.Ny() );
    EXPECT_DOUBLE_EQ( grid.Dx(), q.Dx() );
    EXPECT_ALL_X_EQUAL( grid.getXEdge(lev,i), q.getXEdge(lev,i) );
    EXPECT_ALL_Y_EQUAL( grid.getYEdge(lev,j), q.getYEdge(lev,j) );

    Flux q2( q );
    EXPECT_EQ( nx, q2.Nx() );
    EXPECT_EQ( ny, q2.Ny() );
    EXPECT_DOUBLE_EQ( grid.Dx(), q2.Dx() );
    EXPECT_ALL_X_EQUAL( grid.getXEdge(lev,i), q2.getXEdge(lev,i) );
    EXPECT_ALL_Y_EQUAL( grid.getYEdge(lev,j), q2.getYEdge(lev,j) );
    
    Grid grid2( 2, 3, 1, length, xOffset, yOffset );
    q2.resize( grid2 );
    EXPECT_EQ( 2, q2.Nx() );
    EXPECT_EQ( 3, q2.Ny() );
    EXPECT_DOUBLE_EQ( grid2.Dx(), q2.Dx() );
}

TEST_F(FluxTest, CopyConstructor) {
    Flux h = _f;
    _f = 0.;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), fx(lev,i,j) );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), fy(lev,i,j) );
}

TEST_F(FluxTest, CopyFromDouble) {
    double a = 7;
    _f = a;
    EXPECT_ALL_X_EQUAL( _f(lev,X,i,j), a );
    EXPECT_ALL_Y_EQUAL( _f(lev,Y,i,j), a );
}

TEST_F(FluxTest, UnaryMinus) {
    Flux h = -_f;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), -fx(lev,i,j) );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), -fy(lev,i,j) );
}

TEST_F(FluxTest, PlusEqualsFlux) {
    _g += _f;
    EXPECT_ALL_X_EQUAL( _g(lev,X,i,j), fx(lev,i,j) + gx(lev,i,j) );
    EXPECT_ALL_Y_EQUAL( _g(lev,Y,i,j), fy(lev,i,j) + gy(lev,i,j) );
}

TEST_F(FluxTest, PlusEqualsDouble) {
    _g += 7;
    EXPECT_ALL_X_EQUAL( _g(lev,X,i,j), gx(lev,i,j) + 7 );
    EXPECT_ALL_Y_EQUAL( _g(lev,Y,i,j), gy(lev,i,j) + 7 );
}

TEST_F(FluxTest, MinusEqualsFlux) {
    _g -= _f;
    EXPECT_ALL_X_EQUAL( _g(lev,X,i,j), gx(lev,i,j) - fx(lev,i,j) );
    EXPECT_ALL_Y_EQUAL( _g(lev,Y,i,j), gy(lev,i,j) - fy(lev,i,j) );
}

TEST_F(FluxTest, MinusEqualsDouble) {
    _g -= 7;
    EXPECT_ALL_X_EQUAL( _g(lev,X,i,j), gx(lev,i,j) - 7 );
    EXPECT_ALL_Y_EQUAL( _g(lev,Y,i,j), gy(lev,i,j) - 7 );
}

TEST_F(FluxTest, FluxPlusFlux) {
    Flux h = _f + _g;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), fx(lev,i,j) + gx(lev,i,j) );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), fy(lev,i,j) + gy(lev,i,j) );
}

TEST_F(FluxTest, FluxPlusDouble) {
    Flux h = _f + 7;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), fx(lev,i,j) + 7 );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), fy(lev,i,j) + 7 );        
}

TEST_F(FluxTest, DoublePlusFlux) {
    Flux h = 7 + _f;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), fx(lev,i,j) + 7 );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), fy(lev,i,j) + 7 );        
}

TEST_F(FluxTest, FluxMinusFlux) {
    Flux h = _f - _g;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), fx(lev,i,j) - gx(lev,i,j) );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), fy(lev,i,j) - gy(lev,i,j) );
}

TEST_F(FluxTest, FluxMinusDouble) {
    Flux h = _f - 7;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), fx(lev,i,j) - 7 );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), fy(lev,i,j) - 7 );        
}

TEST_F(FluxTest, DoubleMinusFlux) {
    Flux h = 7 - _f;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), 7 - fx(lev,i,j) );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), 7 - fy(lev,i,j) );        
}

TEST_F(FluxTest, TimesEqualsDouble) {
    _f *= 3;
    EXPECT_ALL_X_EQUAL( _f(lev,X,i,j), 3 * fx(lev,i,j) );
    EXPECT_ALL_Y_EQUAL( _f(lev,Y,i,j), 3 * fy(lev,i,j) );
}

TEST_F(FluxTest, DivEqualsDouble) {
    _f /= 3.;
    EXPECT_ALL_X_EQUAL( _f(lev,X,i,j), fx(lev,i,j) / 3. );
    EXPECT_ALL_Y_EQUAL( _f(lev,Y,i,j), fy(lev,i,j) / 3. );
}

TEST_F(FluxTest, FluxTimesDouble) {
    double a = 5;
    Flux h = _f * a;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), fx(lev,i,j) * a );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), fy(lev,i,j) * a );        
}

TEST_F(FluxTest, DoubleTimesFlux) {
    double a = 5;
    Flux h = a * _f;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), fx(lev,i,j) * a );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), fy(lev,i,j) * a );        
}

TEST_F(FluxTest, FluxDivDouble) {
    double a = 5;
    Flux h = _f / a;
    EXPECT_ALL_X_EQUAL( h(lev,X,i,j), fx(lev,i,j) / a );
    EXPECT_ALL_Y_EQUAL( h(lev,Y,i,j), fy(lev,i,j) / a );        
}   

TEST_F(FluxTest, FluxCoordinates) {
    EXPECT_ALL_X_EQUAL( _f.x(lev,X,i), _grid.getXEdge(lev,i) );
    EXPECT_ALL_X_EQUAL( _f.y(lev,X,j), _grid.getYCenter(lev,j) );
    EXPECT_ALL_Y_EQUAL( _f.x(lev,Y,i), _grid.getXCenter(lev,i) );
    EXPECT_ALL_Y_EQUAL( _f.y(lev,Y,j), _grid.getYEdge(lev,j) );
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
    for (int n=0; n<_ngrid; ++n) {
        for (ind = _f.begin(); ind != _f.end(); ++ind) {
            _f(n,ind) = 35;
        }
    }
    EXPECT_ALL_X_EQUAL( _f(lev,X,i,j), 35 );
    EXPECT_ALL_Y_EQUAL( _f(lev,Y,i,j), 35 );
    for (int lev=0; lev<_ngrid; ++lev) {
        for (ind = _f.begin(); ind != _f.end(); ++ind) {
            EXPECT_DOUBLE_EQ( _f(lev,ind), 35 );
        }
    }
    
    // loop over X alone
    for (int n=0; n<_ngrid; ++n) {
        for (ind = _f.begin(X); ind != _f.end(X); ++ind) {
            _f(n,ind) = 53;
        }
    }
    EXPECT_ALL_X_EQUAL( _f(lev,X,i,j), 53 );
    EXPECT_ALL_Y_EQUAL( _f(lev,Y,i,j), 35 );

    // loop over Y alone
    for (int n=0; n<_ngrid; ++n) {
        for (ind = _f.begin(Y); ind != _f.end(Y); ++ind) {
            _f(n,ind) = 350;
        }
    }
    EXPECT_ALL_X_EQUAL( _f(lev,X,i,j), 53 );
    EXPECT_ALL_Y_EQUAL( _f(lev,Y,i,j), 350 );

}

TEST_F(FluxTest, GetIndex) {
    EXPECT_ALL_X_EQUAL( _f(lev,_f.getIndex(X,i,j)), _f(lev,X,i,j) );
    EXPECT_ALL_Y_EQUAL( _f(lev,_f.getIndex(Y,i,j)), _f(lev,Y,i,j) );
}

TEST_F(FluxTest, GridValues) {
    EXPECT_ALL_X_EQUAL( _f.x(lev,_f.getIndex(X,i,j)), _f.x(lev,X,i) );
    EXPECT_ALL_Y_EQUAL( _f.x(lev,_f.getIndex(Y,i,j)), _f.x(lev,Y,i) );    
    EXPECT_ALL_X_EQUAL( _f.y(lev,_f.getIndex(X,i,j)), _f.y(lev,X,j) );
    EXPECT_ALL_Y_EQUAL( _f.y(lev,_f.getIndex(Y,i,j)), _f.y(lev,Y,j) );    
}

TEST_F(FluxTest, UniformFlow) {
    double mag = 3.;
    double angle = 0.;
    Flux q = Flux::UniformFlow( _grid, mag, angle );
    EXPECT_ALL_X_EQUAL( q(lev,X,i,j), mag * _grid.Dx(lev) );
    EXPECT_ALL_Y_EQUAL( q(lev,Y,i,j), 0. );

    double pi = 4. * atan(1.);
    angle = pi/2.;
    q = Flux::UniformFlow( _grid, mag, angle );
    EXPECT_ALL_X_NEAR( q(lev,X,i,j), 0., 1e-15 );
    EXPECT_ALL_Y_EQUAL( q(lev,Y,i,j), mag * _grid.Dx(lev) );
}

#undef EXPECT_ALL_X_EQUAL
#undef EXPECT_ALL_Y_EQUAL
