#include "BoundaryVector.h"
#include <gtest/gtest.h>

using namespace ibpm;

#define EXPECT_ALL_EQUAL(a,b) \
    for ( int i=0; i<_n; ++i ) { \
        EXPECT_DOUBLE_EQ( (a), (b) ); \
    }

class BoundaryVectorTest : public testing::Test {
protected:
    BoundaryVectorTest() :
        _n(10),
        _f(_n),
        _g(_n)
    {
        for (int i=0; i<_n; ++i) {
            _f(X,i) = fx(i);
            _f(Y,i) = fy(i);
            _g(X,i) = gx(i);
            _g(Y,i) = gy(i);
        }
    }

    double fx(int i) { return 2*i; }
    
    double fy(int i) { return 3*i; }

    double gx(int i) { return _n-i; }

    double gy(int i) { return 2*_n-i; }

    int _n;
    BoundaryVector _f;
    BoundaryVector _g;
};

TEST_F(BoundaryVectorTest, TestAssignment) {
    int i = _n/2;
    _f(X,i) = 73;
    EXPECT_DOUBLE_EQ(_f(X,i), 73 );
    _f(Y,i) = 23;
    EXPECT_DOUBLE_EQ(_f(Y,i), 23 );
}

TEST_F(BoundaryVectorTest, testCopyConstructor) {
    BoundaryVector h = _f;
    _f = 0.;
    EXPECT_ALL_EQUAL( h(X,i), fx(i) );
    EXPECT_ALL_EQUAL( h(Y,i), fy(i) );
}

TEST_F(BoundaryVectorTest, testCopyFromDouble) {
    double a = 7;
    _f = a;
    EXPECT_ALL_EQUAL( _f(X,i), a );
    EXPECT_ALL_EQUAL( _f(Y,i), a );
}

TEST_F(BoundaryVectorTest, testUnaryMinus) {
    BoundaryVector h = -_f;
    EXPECT_ALL_EQUAL( h(X,i), -fx(i) );
    EXPECT_ALL_EQUAL( h(Y,i), -fy(i) );
}

TEST_F(BoundaryVectorTest, PlusEqualsBoundaryVector) {
    _g += _f;
    EXPECT_ALL_EQUAL( _g(X,i), fx(i) + gx(i) );
    EXPECT_ALL_EQUAL( _g(Y,i), fy(i) + gy(i) );
}

TEST_F(BoundaryVectorTest, MinusEqualsBoundaryVector) {
    _g -= _f;
    EXPECT_ALL_EQUAL( _g(X,i), gx(i) - fx(i) );
    EXPECT_ALL_EQUAL( _g(Y,i), gy(i) - fy(i) );
}

TEST_F(BoundaryVectorTest, BoundaryVectorPlusBoundaryVector) {
    BoundaryVector h = _f + _g;
    EXPECT_ALL_EQUAL( h(X,i), fx(i) + gx(i) );
    EXPECT_ALL_EQUAL( h(Y,i), fy(i) + gy(i) );
}

TEST_F(BoundaryVectorTest, BoundaryVectorMinusBoundaryVector) {
    BoundaryVector h = _f - _g;
    EXPECT_ALL_EQUAL( h(X,i), fx(i) - gx(i) );
    EXPECT_ALL_EQUAL( h(Y,i), fy(i) - gy(i) );
}

TEST_F(BoundaryVectorTest, TimesEqualsDouble) {
    _f *= 3;
    EXPECT_ALL_EQUAL( _f(X,i), 3 * fx(i) );
    EXPECT_ALL_EQUAL( _f(Y,i), 3 * fy(i) );
}

TEST_F(BoundaryVectorTest, DivEqualsDouble) {
    _f /= 3.;
    EXPECT_ALL_EQUAL( _f(X,i), fx(i) / 3. );
    EXPECT_ALL_EQUAL( _f(Y,i), fy(i) / 3. );
}

TEST_F(BoundaryVectorTest, TimesDouble) {
    double a = 7;
    BoundaryVector h = _f * a;
    EXPECT_ALL_EQUAL( h(X,i), fx(i) * a );
    EXPECT_ALL_EQUAL( h(Y,i), fy(i) * a );
}

TEST_F(BoundaryVectorTest, DoubleTimes) {
    double a = 7;
    BoundaryVector h = a * _f;
    EXPECT_ALL_EQUAL( h(X,i), fx(i) * a );
    EXPECT_ALL_EQUAL( h(Y,i), fy(i) * a );
}

TEST_F(BoundaryVectorTest, DivDouble) {
    double a = 7;
    BoundaryVector h = _f / a;
    EXPECT_ALL_EQUAL( h(X,i), fx(i) / a );
    EXPECT_ALL_EQUAL( h(Y,i), fy(i) / a );
}

TEST_F(BoundaryVectorTest, IndexCount) {
    BoundaryVector::index ind;
    int count=0;
    for (ind = _f.begin(); ind != _f.end(); ++ind) {
        ++count;
    }
    EXPECT_EQ(count, _n*2);

    count=0;
    for (ind = _f.begin(X); ind != _f.end(X); ++ind) {
        ++count;
    }
    EXPECT_EQ(count, _n);
    
    count=0;
    for (ind = _f.begin(Y); ind != _f.end(Y); ++ind) {
        ++count;
    }
    EXPECT_EQ(count, _n);
}

TEST_F(BoundaryVectorTest, IndexReference) {
    BoundaryVector::index ind;

    // loop over all values
    for (ind = _f.begin(); ind != _f.end(); ++ind) {
        _f(ind) = 35;
    }
    EXPECT_ALL_EQUAL( _f(X,i), 35 );
    EXPECT_ALL_EQUAL( _f(Y,i), 35 );
    for (ind = _f.begin(); ind != _f.end(); ++ind) {
        EXPECT_DOUBLE_EQ( _f(ind), 35 );
    }
    
    // loop over X alone
    for (ind = _f.begin(X); ind != _f.end(X); ++ind) {
        _f(ind) = 53;
    }
    EXPECT_ALL_EQUAL( _f(X,i), 53 );
    EXPECT_ALL_EQUAL( _f(Y,i), 35 );

    // loop over Y alone
    for (ind = _f.begin(Y); ind != _f.end(Y); ++ind) {
        _f(ind) = 350;
    }
    EXPECT_ALL_EQUAL( _f(X,i), 53 );
    EXPECT_ALL_EQUAL( _f(Y,i), 350 );

}

TEST_F(BoundaryVectorTest, GetIndex) {
    BoundaryVector::index ind;
    for (int i = 0; i < _n; ++i) {
        ind = _f.getIndex(X, i);
        _f(ind) = gx(i);
        ind = _f.getIndex(Y, i);
        _f(ind) = gy(i);
    }
    EXPECT_ALL_EQUAL( _f(X,i), gx(i) );
    EXPECT_ALL_EQUAL( _f(Y,i), gy(i) );    
}

TEST_F(BoundaryVectorTest, NoPoints) {
    BoundaryVector b(0);
    EXPECT_EQ( b.getNumPoints(), 0 );
}

#undef EXPECT_ALL_EQUAL
