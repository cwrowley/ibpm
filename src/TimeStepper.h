#ifndef _TIMESTEPPER_H_
#define _TIMESTEPPER_H_

#include <string>
using std::string;

namespace ibpm {

class Model;
class ProjectionSolver;
class State;

/*!
    \file TimeStepper.h
    \class TimeStepper

    \brief Abstract base class for advancing a flow field forward in time
    
    The governing equations are in the form
    \f{align}
       \frac{d\gamma}{dt} + Bf &= \alpha L\gamma + N(q)\\
       C\gamma &= b
    \f}
    where \f$ \gamma \f$ is a Scalar and \f$ f \f$ and \f$ b \f$ are
    BoundaryVectors.  \f$ L \f$ is the discrete Laplacian.  The operators
    B, C, and N, as well as the constant \f$ \alpha \f$ and the constraint
    values \f$ b \f$, are determined by an associated Model.
    
    \author Clancy Rowley
    \author $LastChangedBy$
    \date  2 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class TimeStepper {

//
// Public methods
//
public:
    
    /*! \brief Setup routines necessary to step the solution forward in time.
    
    Note that creation of the ProjectionSolver should be deferred to the
    subclasses, but determination of which type of solver to instantiate
    is handled by the base class.
    \param[in] name     The name of the solver (e.g. Explicit Euler)
    \param[in] model    The associated Model instance
    \param[in] timestep Timestep to use for the advance() routine
    */
    TimeStepper(string name, Grid& grid, Model& model, double timestep);

    virtual ~TimeStepper() {}

    /*! \brief Perform initialization, if necessary
    
    Calls init() on the associated Model, in addition to
    initializations needed within TimeStepper
    */
    virtual void init();
    
    /*! \brief Load the state of the solver from a file (see saveState)

    Return true if successful
    */
    virtual bool load(const string& filename);
    
    /*! \brief Save the state of the solver to the specified file

    For instance, this might be a Cholesky factorization saved to avoid
    recomputing, or information at a previous timestep for AdamsBashforth.
    Can be used in place of init() (if successful)
    Return true if successful
    */
    virtual bool save(const string& filename);

    /*! \brief Return the name of the timestepping scheme, as a string
    */
    string getName();

    /*! \brief Advance the state forward in time.
    
    Pure virtual method: must be overridden by subclasses
    \param[inout] x State vector to be stepped forward in time
    */
    virtual void advance(State& x) = 0;

//
// Protected data and methods
//
protected:

    const string _name;
    Grid& _grid;
    Model& _model;
    double _timestep;

    /*! \brief Return a new ProjectionSolver of the appropriate type
    
    \param[in] beta     Parameter passed to the appropriate ProjectionSolver
                        (usually the timestep)
    */
    ProjectionSolver* createSolver(double beta);

};

}

#endif /* _TIMESTEPPER_H_ */
