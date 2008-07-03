#ifndef _PROJECTIONSOLVER_H_
#define _PROJECTIONSOLVER_H_

/**
 * \file ProjectionSolver.h
 * \class ProjectionSolver
 *
 * \brief Solve systems arising in a fractional step method.
 * 
 * Solve a system of the form
 * \f[
 *    (1 - \frac{\alpha}{2}L)x + \alpha Bf = a
 * \f]\f[
 *    Cx = b
 * \f]
 * for f and x, using the following algorithm (a fractional step method):
 * 
 * \f{align*}{
 *    Ax^* &= a\\
 *    C A^{-1}Bf &= Cx^* - b\\
 *    x &= x^* - A^{-1}Bf
 * \f}
 * where \f$ A = 1-(\alpha/2)L \f$.
 * 
 * Note that this form of the equations arises for all of the timesteppers
 * considered (so far).  The parameter \f$ \alpha \f$ is fixed, and the
 * operators L, B, and C are determined by an associated NavierStokesModel.
 * 
 * This is an abstract base class, and may not be instantiated directly.
 * 
 * \author Clancy Rowley
 * \date  3 Jul 2008
 *
 * $Revision: $
 * $LastChangedDate: $
 * $LastChangedBy: $
 * $HeadURL: $
 */

class ProjectionSolver {
  public:
	ProjectionSolver(const NavierStokesModel& model, double alpha) :
		_alpha(alpha), _model(&model);
	~ProjectionSolver();
	void solve(const Scalar& a, const BoundaryVector& b, Scalar& gamma, BoundaryVector& f);
  protected:
	void Ainv(const Scalar& x, Scalar& y);
	void B(const BoundaryVector& f, Scalar& y);
	void C(const Scalar& x, BoundaryVector& y);
	void M(const BoundaryVector& f, BoundaryVector& y);
	virtual void Minv(const BoundaryVector& b, BoundaryVector& f) = 0;
	Geometry* getGeometry();
  private:
	double _alpha;
	NavierStokesModel* _model;
	
}

#endif /* _PROJECTIONSOLVER_H_ */
