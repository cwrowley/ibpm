#ifndef _CHOLESKYSOLVER_H_
#define _CHOLESKYSOLVER_H_

#include "ProjectionSolver.h"
#include "Array.h"
#include <string>

using Array::array1;
using Array::array2;

namespace ibpm {

/*!
    \file CholeskySolver.h
    \class CholeskySolver

    \brief Subclass of ProjectionSolver, in which the system M f = b is solved 
    directly, using Cholesky factorization.
 
    Here M is the matrix 
        \beta C A^{-1} B 
    where
        A = (1 + \alpha\beta / 2) * Laplacian
 
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
        const Grid& grid,
        const NavierStokesModel& model,
        double beta
    );
    
    ~CholeskySolver();

    /// \brief Compute the Cholesky decomposition of M
    void init();

    /// \brief Load a Cholesky decomposition from the specified file.
    /// Returns true if successful
    bool load(const std::string& filename);
    
    /// \brief Save a Cholesky decomposition to the specified file, 
    /// overwriting if necessary.
    /// Returns true if successful
    bool save(const std::string& filename);
    
protected:
    /// \brief Solve Mf = b for f, using the Cholesky factorization of M. 
    /// Assumes M is symmetric.
    void Minv(
        const BoundaryVector& b,
        BoundaryVector& x
    );

private:
    int _numPoints;  // number of points in the geometry
    int _size;       // size of the vectors: numPoints * 2
    const double _alphaBeta;    // keep a local copy of alpha * beta, as a
                                // check when reading in Cholesky factorizations
                                // from files
    array2<double> _lower;
    array1<double> _diagonal;
    void computeMatrixM( array2<double>& M );
    void computeFactorization( const array2<double>& M );
    bool _hasBeenInitialized;
};

} // namespace ibpm

#endif /* _CHOLESKYSOLVER_H_ */
