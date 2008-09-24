#ifndef _MODEL_H_
#define _MODEL_H_

#include "Scalar.h"
#include "BoundaryVector.h"
#include "State.h"

namespace ibpm {
    

/*!
    \file Model.h
    \class Model

    \brief Abstract base class defining a Model for use with TimeStepper
    and ProjectionSolver
    
    Subclasses define operators for models of the form
    \f{align}
       \frac{d\gamma}{dt} + Bf &= \alpha L\gamma + N(q)\\
       C\gamma &= b
    \f}
    where \f$ \gamma \f$ is a Scalar and \f$ f \f$ and \f$ b \f$ are
    BoundaryVectors, and \f$ q \f$ is a State.
    \f$ L \f$ is the discrete Laplacian, and Model subclasses are 
    responsible for determining operators B, C, and N, as well as the
    constant \f$ \alpha \f$.

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 18 Sep 2008
    \date $LastChangedDate$
    \version $Revision$
*/

    class Model {
    public:
        Model() {}
        virtual ~Model() {}
        virtual void init() {}
        
        /// \brief Return true if the operators depend on time (e.g. moving body)
        virtual bool isTimeDependent() const = 0;
        
        /// \brief Return the number of points in constraint BoundaryVector b
        virtual int getNumPoints() const = 0;

        /// \brief Rreturn the right-hand side b of the constraint equation
        virtual BoundaryVector getConstraints() const = 0;
        
        /// \brief update operators for time-dependent models
        virtual void updateOperators( double time ) = 0;

        virtual void B( const BoundaryVector& f, Scalar& gamma ) const = 0;
        virtual void C( const Scalar& gamma, BoundaryVector& f ) const = 0;
        virtual Scalar N( const State& x ) const = 0;
        virtual double getAlpha() const = 0;
        
        // \brief Compute any quantities in the state vector that depend on
        // other state variables and the model (e.g. for Navier-Stokes,
        // compute the flux from the circulation)
        virtual void refreshState( State& x ) const {}
    };

} // namespace ibpm

#endif /* _MODEL_H_ */
