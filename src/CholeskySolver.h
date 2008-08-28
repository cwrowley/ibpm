#ifndef _CHOLESKYSOLVER_H_
#define _CHOLESKYSOLVER_H_

#include "ProjectionSolver.h"

/*!
    \file CholeskySolver.h
    \class CholeskySolver

    \brief Subclass of ProjectionSolver, in which the system M f = b is solved 
    directly, using Cholesky factorization.

    This class assumes the matrix M is symmetric.  For now, it computes the 
    Cholesky factorization when it is instantiated, but this should be changed 
    to do this explicitly, or load a previously computed factorization from a 
    file.

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 28 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class CholeskySolver : public ProjectionSolver {
public:
    
    CholeskySolver(
        const NavierStokesModel& model,
        double alpha
    );
    
    ~CholeskySolver();

    /// \brief Compute the Cholesky decomposition of M
    void computeCholesky();

    /// \brief Load a Cholesky decomposition from the specified file.
    /// Returns true if successful
    bool loadCholesky(const string& filename);
    
    /// \brief Save a Cholesky decomposition to the specified file, 
    /// overwriting if necessary.
    /// Returns true if successful
    bool saveCholesky(const string& filename);
    
protected:
    /// \brief Solve Mf = b for f, using the Cholesky factorization of M. 
    /// Assumes M is symmetric.
    void Minv(
        const BoundaryVector& b,
        BoundaryVector& x
    );

private:
    double* _M;
};


#endif /* _CHOLESKYSOLVER_H_ */
