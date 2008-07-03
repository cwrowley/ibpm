#ifndef _PROJECTIONSOLVER_H_
#define _PROJECTIONSOLVER_H_

/*!
\file ProjectionSolver.h
\class ProjectionSolver

\brief Solve systems arising in a fractional step method.

Solve a system of the form
\f{align}
   (1 - \frac{\alpha}{2}L)x + \alpha Bf &= a\\
   Cx &= b
\f}
for f and x, using the following algorithm (a fractional step method):

\f{align}
   Ax^* &= a\\
   C A^{-1}Bf &= Cx^* - b\\
   x &= x^* - A^{-1}Bf
\f}

where
\f[
   A = 1-\frac{\alpha}{2}L.
\f]

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

ProjectionSolver(const NavierStokesModel& model, double alpha) :
		_alpha(alpha), _model(&model);
~ProjectionSolver();

/*! \brief Solve for \a gamma and \a f using a fractional step method.
Solves equations (1-2) using the algorithm (3-5).
\param[in] a	Right-hand side of equation (1)
\param[in] b	Right-hand side of constraint equation (2)
\param[out] gamma	Circulation after solution
\param[out] f	Boundary force after solution
*/
void solve(const Scalar& a, const BoundaryVector& b, Scalar& gamma, BoundaryVector& f);

//
// Protected methods
//
protected:

/// Solve \f$ y = A^{-1} b \f$.
void Ainv(const Scalar& x, Scalar& y);

/// Compute \a y = B(\a f)
void B(const BoundaryVector& f, Scalar& y);

/// Compute \a y = C(\a x)
void C(const Scalar& x, BoundaryVector& y);

/// Compute \a y = M(\a f), where \f$ M = C A^{-1} B \f$.
void M(const BoundaryVector& f, BoundaryVector& y);

/// Compute \f$ f = M^{-1} b \f$.
virtual void Minv(const BoundaryVector& b, BoundaryVector& f) = 0;

Geometry* getGeometry();

//
// Private data
//
private:
	double _alpha;
	NavierStokesModel* _model;
	
}

#endif /* _PROJECTIONSOLVER_H_ */
