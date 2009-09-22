#ifndef _PROJECTIONSOLVER_H_
#define _PROJECTIONSOLVER_H_

#include "Scalar.h"
#include "BoundaryVector.h"
#include "NavierStokesModel.h"
#include "EllipticSolver.h"
#include <string>
using std::string;

namespace ibpm {

/*!
\file ProjectionSolver.h
\class ProjectionSolver

\brief Solve systems arising in a fractional step method.

Solve a system of the form
\f{align}
   (1 - \frac{\alpha\beta}{2}L )x + \beta B f &= a\\
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
   A &= 1-\frac{\alpha\beta}{2}L\\
   B_1 &= \beta B
\f}

Note that this form of the equations arises for all of the timesteppers
considered (so far).  The parameter \f$ \beta \f$ is fixed, L is the discrete
Laplacian, and the operators B, and C and the parameter \f$ \alpha \f$ are
determined by an associated Model.

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
    ProjectionSolver(const Grid& grid, const NavierStokesModel& model, double beta);

    /// Destructor
    virtual ~ProjectionSolver();

    /// Perform initialization, if needed
    /// NOTE:  Does not call init() on the Model model
    virtual void init();

    /// \brief Save information to a file, if there is any to write.
    /// Return true if successful
    virtual bool save(const string& basename);
    
    /// \brief Load information from a file
    /// Can be used in place of init() (if successful)
    /// Return true if successful
    virtual bool load(const string& basename);
    
    /*! \brief Solve for \a omega and \a f using a fractional step method.
    Solves equations (1-2) using the algorithm (3-5).
    
    Assumes the Model has been initialized with init()
    
    \param[in] a	Right-hand side of equation (1)
    \param[in] b	Right-hand side of constraint equation (2)
    \param[out] omega	Vorticity after solution
    \param[out] f	Boundary force after solution
    */
    void solve(
        const Scalar& a,
        const BoundaryVector& b,
        Scalar& omega,
        BoundaryVector& f
    );

//
// Protected methods
//
protected:

    /// Solve \f$ y = A^{-1} b \f$.
    void Ainv( const Scalar& b, Scalar& y );

    /// Compute \a y = B(\a f)
    void B( const BoundaryVector& f, Scalar& y );
    
    /// Compute \a y = C(\a x)
    void C( const Scalar& x, BoundaryVector& y );
    
    /// Compute \a y = M(\a f), where \f$ M = C A^{-1} B \f$.
    void M( const BoundaryVector& f, BoundaryVector& y );
    BoundaryVector M( const BoundaryVector& f );
    
    /// Compute \f$ x = M^{-1} b \f$.
    virtual void Minv( const BoundaryVector& b, BoundaryVector& x ) = 0;
    
//
// Private data
//
private:
	double _beta;
    const Grid _grid;
	const NavierStokesModel& _model;
    HelmholtzSolver _helmholtz;
};

} // namespace ibpm

#endif /* _PROJECTIONSOLVER_H_ */
