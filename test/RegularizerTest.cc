#include "Regularizer.h"
#include "Grid.h"
#include "Scalar.h"
#include "RigidBody.h"
#include "Geometry.h"
#include "VectorOperations.h"
#include <gtest/gtest.h>
#include <iostream>

using namespace std;
using namespace ibpm;

namespace {

double tol = 1e-14;
    
class RegularizerTest : public testing::Test {
protected:
    RegularizerTest() : 
        _nx(10),
        _ny(20),
        _ngrid(1),
        _grid(_nx, _ny, _ngrid, 2, -1, -2),
        _dx( _grid.Dx() ),
        _u1(1),
        _u2(_grid)
    {
        // Create a geometry with one point
        RigidBody body;
        body.addPoint(0,0);
        _geom.addBody(body);
        _regularizer = new Regularizer(_grid, _geom);
        _regularizer->update();

        _u1(X,0) = 3;
        _u1(Y,0) = 5;
        _u2 = _regularizer->toFlux(_u1);        
    }

    // data
    int _nx;
    int _ny;
    int _ngrid;
    Grid _grid;
    double _dx;
    BoundaryVector _u1;
    Flux _u2;
    Geometry _geom;
    Regularizer* _regularizer;
};

// Properties of phi, from Roma 1999 (Candidates for tests):
//
//     1. phi(r) is continuous
//         Derivative:
//             -r/(\sqrt(-3r^2 + 1)), |r| < 0.5
//             -1/2 * r/abs(r) - 1/( 2 * sqrt(-3(1-abs(r))^2 + 1) ) *
//                  6(1-abs(r)) * r/abs(r), 0.5 < |r| < 1.5;
//         max slope is 1
//         could test that abs(delta(x+eps) - delta(x)) <= eps  for all x
//         BUT, should only test the public interface, not the delta function
//
//     2. phi(r) = 0 for |r| >= 1.5
//         - test delta function directly
//         - smoothing: flux is zero for points > 1.5 pts away
//         - interpolating: boundary points are zero if > 1.5 pts away from
//           nonzero fluxes
//
//     3. \sum phi(r-i) = 1 for all r
//         - Flux values should add up to sum of boundary points 
//
//     4. \sum (r-i) phi(r-i) = 0 for all r
//         - Moment of a single point should be zero
//
//     5. phi(r-i)^2 = 1/2 for all r
//         - (boundary -> grid -> boundary) should give 1/4 * original value

TEST_F(RegularizerTest, FluxSumsToOne) {
    Scalar one( _grid );
    Scalar zero( _grid );
    Flux q( _grid );

    one = 1.;
    zero = 0.;
    
    // Sum X component
    VelocityToFlux( one, zero, q );
    double sum = InnerProduct( q, _u2 );
    // Normalize by dx^2: inner product is
    // \int (u1 * u2) dx dy \approx \sum (u1 * u2) * dx^2
    sum /= _dx * _dx;
    EXPECT_DOUBLE_EQ( _u1(X,0), sum );

    // Sum Y component
    VelocityToFlux( zero, one, q );
    sum = InnerProduct( q, _u2);
    // Normalize by dx^2, as before
    sum /= _dx * _dx;
    EXPECT_DOUBLE_EQ( _u1(Y,0), sum );
}

TEST_F(RegularizerTest, MomentIsZero) {
    Flux::index ind;
    double moment=0;
    for (ind = _u2.begin(X); ind != _u2.end(X); ++ind) {
        moment += _u2.x(0,ind) * _u2(0,ind);
    }
    EXPECT_NEAR(moment, 0, tol);
    
    moment = 0;
    for (ind = _u2.begin(X); ind != _u2.end(X); ++ind) {
        moment += _u2.x(0,ind) * _u2(0,ind);
    }
    EXPECT_NEAR(moment, 0, tol);
}

TEST_F(RegularizerTest, CountNumberNonzeroPoints) {
    Flux::index ind;
    int count=0;
    for (ind = _u2.begin(X); ind != _u2.end(X); ++ind) {
        if (fabs(_u2(0,ind)) > tol) {
#ifdef DEBUG
            cout << "(" << _u2.x(ind) << "," << _u2.y(ind) << ") : " <<
                _u2(ind) << endl;
#endif
            ++ count;
        }
    }
    // Should have 6 nonzero points:
    //   x = {-0.2, 0, 0.2}
    //   y = {-0.1, 0.1}
    EXPECT_EQ(count, 6);

    count=0;
    for (ind = _u2.begin(Y); ind != _u2.end(Y); ++ind) {
        if (fabs(_u2(0,ind)) > tol) {
#ifdef DEBUG
            cout << "(" << _u2.x(ind) << "," << _u2.y(ind) << ") : " << 
                _u2(ind) << endl;
#endif
            ++ count;
        }
    }
    // Should have 6 nonzero points:
    //   x = {-0.1, 0.1}
    //   y = {-0.2, 0, 0.2}
    EXPECT_EQ(count, 6);
    
}

TEST_F(RegularizerTest, FluxZeroAwayFromBoundary) {
    Flux::index ind;
    double supportDistance = 1.5 * _grid.Dx();
    
    for (ind = _u2.begin(); ind != _u2.end(); ++ind) {
        // if point is outside of region of support of the origin
        if (fabs(_u2.x(0,ind)) > supportDistance || 
            fabs(_u2.y(0,ind)) > supportDistance) {
            EXPECT_DOUBLE_EQ(_u2(0,ind), 0);
        }
        if (fabs(_u2.x(0,ind)) + tol < supportDistance && 
            fabs(_u2.y(0,ind)) + tol < supportDistance) {
#ifdef DEBUG
            cout << "(" << _u2.x(ind) << "," << _u2.y(ind) << ") : " <<
               _u2(ind) << endl;
#endif
            EXPECT_GT(fabs(_u2(0,ind)), 0);
        }
    }
}

TEST_F(RegularizerTest, InterpolateConstant) {
    Scalar u(_grid);
    Scalar v(_grid);
    Flux   q(_grid);
    u = 7;
    v = 9;
    VelocityToFlux( u, v, q );
    BoundaryVector f = _regularizer->toBoundary( q );
    
    EXPECT_NEAR(f(X,0), 7, tol);
    EXPECT_NEAR(f(Y,0), 9, tol);
}

TEST_F(RegularizerTest, SmoothThenInterpolateEqualsQuarter) {
    BoundaryVector f = _regularizer->toBoundary(_u2);
    
    EXPECT_NEAR( f(X,0), _u1(X,0) * 0.25, tol );
    EXPECT_NEAR( f(Y,0), _u1(Y,0) * 0.25, tol );
}

} // namespace
