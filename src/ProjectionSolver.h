#ifndef _PROJECTIONSOLVER_H_
#define _PROJECTIONSOLVER_H_

#include "Scalar.h"
#include "BoundaryVector.h"
#include "NavierStokesModel.h"

namespace ibpm {

/*!
\file ProjectionSolver.h
\class ProjectionSolver

\brief Solve systems arising in a fractional step method.

Solve a system of the form
\f{align}
   (1 - \frac{\alpha}{2}L )x + \alpha B f &= a\\
   C x &= b
\f}
for \a f and \a x, using the following algorithm (a fractional step method):

\f{align}
   Ax^* &= a\\
   C A^{-1}B_1 f &= Cx^* - b\\
   x &= x^* - A^{-1}B_1 f
\f}

where
\f{align}
   A &= 1-\frac{\alpha}{2}L\\
   B_1 &= \alpha B
\f}

Note that this form of the equations arises for all of the timesteppers
considered (so far).  The parameter \f$ \alpha \f$ is fixed, and the
operators L, B, and C are determined by an associated NavierStokesModel.

This is an abstract base class, and may not be instantiated directly.

\author Clancy Rowley
\date  3 Jul 2008

$Revision$
$LastChangedDate$
$LastChangedBy$
$HeadURL$
*/

class ProjectionSolver {

//
// Public methods
//
public:

    /// Constructor
    ProjectionSolver(const NavierStokesModel& model, double alpha);

    /// Destructor
    virtual ~ProjectionSolver();

    /*! \brief Solve for \a gamma and \a f using a fractional step method.
    Solves equations (1-2) using the algorithm (3-5).
    \param[in] a	Right-hand side of equation (1)
    \param[in] b	Right-hand side of constraint equation (2)
    \param[out] gamma	Circulation after solution
    \param[out] f	Boundary force after solution
    */

    void solve(
        const Scalar& a,
        const BoundaryVector& b,
        Scalar& gamma,
        BoundaryVector& f
    );

//
// Protected methods
//
protected:

    /// Solve \f$ y = A^{-1} b \f$.
    inline Scalar Ainv(const Scalar& b) {
        Scalar Sb = _model->S( b );
        Sb *= _eigenvaluesOfAinv;
        Scalar y = _model->Sinv( Sb );
        return y;
    }
    
    /// Solve \f$ y = A^{-1} b \f$.
    inline Scalar Ainv(const Scalar& b, Scalar& y ) {
        y = _model->S( b );
        y *= _eigenvaluesOfAinv;
        y = _model->Sinv( y );
        return y;
    }

    /// Compute \a y = B(\a f)
    inline Scalar B(const BoundaryVector& f) {
        Scalar y = _model->B( f );
        y *= _alpha;
        return y;
    }
    
    /// Compute \a y = C(\a x)
    inline BoundaryVector C(const Scalar& x) {
        BoundaryVector y = _model->C( x );
        return y;
    }
    
    /// Compute \a y = C(\a x)
    inline void C(const Scalar& x, BoundaryVector& y ) {
        y = _model->C( x );
    }
    
    /// Compute \a y = M(\a f), where \f$ M = C A^{-1} B \f$.
    inline BoundaryVector M(const BoundaryVector& f) {
        Scalar Bf = B(f);
        Scalar Ainv_B_f = Ainv(Bf);
        BoundaryVector y = C(Ainv_B_f);
        return y;
    }
    
    /// Compute \a y = M(\a f), where \f$ M = C A^{-1} B \f$.
    inline void M(const BoundaryVector& f, BoundaryVector& y ) {
        Scalar Bf = B( f );
        Ainv( Bf, Bf );
        C( Bf, y );
    }    

    /// Compute \f$ x = M^{-1} b \f$.
    virtual void Minv(
        const BoundaryVector& b,
        BoundaryVector& x
    ) = 0;
    
    /// Return a pointer to the associated geometry
    inline const Geometry* getGeometry() {
        return _model->getGeometry();
    }

//
// Private data
//
private:
	double _alpha;
	const NavierStokesModel* _model;
    Scalar _eigenvaluesOfAinv;	
};

} // namespace ibpm

#endif /* _PROJECTIONSOLVER_H_ */
