/// \file exceptions.h
/// \brief Custom exceptions for the \c Newton_GMRES suite.
/// \author Kevin K. Chen (Princeton University)

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>
#include <string>

// ************************************************************************** //

namespace arnoldi_ns {

using std::runtime_error;

/// Exception thrown when something goes wrong in \c ArnoldiDecomp.
class ArnoldiError : public runtime_error {
public:
    /// Construct with a message.
    explicit ArnoldiError(const std::string& msg) : runtime_error{msg} {}
    ArnoldiError(const ArnoldiError&) = default; ///< Copy constructor.
    ArnoldiError(ArnoldiError&&) = default; ///< Move constructor.

    ~ArnoldiError() noexcept override {} ///< Destructor (empty).

private:
    void operator=(const ArnoldiError&) = delete;
    void operator=(ArnoldiError&&) = delete;
};

} // End namespace arnoldi_ns

// ************************************************************************** //

namespace continuation_ns {

using std::runtime_error;

/// Exception thrown when something goes wrong in \c Continuation.
class ContinuationError : public runtime_error {
public:
    /// Construct with a message.
    explicit ContinuationError(const std::string& msg) : runtime_error{msg} {}
    /// Copy constructor.
    ContinuationError(const ContinuationError&) = default;
    ContinuationError(ContinuationError&&) = default; ///< Move constructor.

    ~ContinuationError() noexcept override {} ///< Destructor (empty).

private:
    void operator=(const ContinuationError&) = delete;
    void operator=(ContinuationError&&) = delete;
};

} // End namespace continuation_ns

// ************************************************************************** //

/// The namespace for auxiliary constructs.
namespace aux_ns {

using std::runtime_error;

/// Exception thrown when the user gives an illegal value for a parameter.
class InvalidValue : public runtime_error {
public:
    /// Construct with a message.
    explicit InvalidValue(const std::string& msg) : runtime_error{msg} {}
    InvalidValue(const InvalidValue&) = default; ///< Copy constructor.
    InvalidValue(InvalidValue&&) = default; ///< Move constructor.

    ~InvalidValue() noexcept override {} ///< Destructor (empty).

private:
    void operator=(const InvalidValue&) = delete;
    void operator=(InvalidValue&&) = delete;
};

} // End namespace continuation_ns

// ************************************************************************** //

namespace linear_solver_ns {

using std::runtime_error;

/// Exception thrown when something goes wrong in \c LinearSolver objects.
class LinearSolverError : public runtime_error {
public:
    /// Construct with a message.
    explicit LinearSolverError(const std::string& msg) : runtime_error{msg} {}
    /// Copy constructor.
    LinearSolverError(const LinearSolverError&) = default;
    LinearSolverError(LinearSolverError&&) = default; ///< Move constructor.

    ~LinearSolverError() noexcept override {} ///< Destructor (empty).

private:
    void operator=(const LinearSolverError&) = delete;
    void operator=(LinearSolverError&&) = delete;
};

} // End namespace continuation_ns

// ************************************************************************** //

namespace newton_ns {

using std::runtime_error;

/// Exception thrown when something goes wrong in \c Newton objects.
class NewtonError : public runtime_error {
public:
    /// Construct with a message.
    explicit NewtonError(const std::string& msg) : runtime_error{msg} {}
    NewtonError(const NewtonError&) = default; ///< Copy constructor.
    NewtonError(NewtonError&&) = default; ///< Copy constructor.

    ~NewtonError() noexcept override {} ///< Destructor (empty).

private:
    void operator=(const NewtonError&) = delete;
    void operator=(NewtonError&&) = delete;
};

} // End namespace newton_ns

// ************************************************************************** //

namespace aux_ns {

using std::logic_error;

/// Exception thrown when attempting to use an order that is not implemented.
class OrderError : public logic_error {
public:
    /// Construct with a message.
    explicit OrderError(const std::string& msg) : logic_error{msg} {}
    OrderError(const OrderError&) = default; ///< Copy constructor.
    OrderError(OrderError&&) = default; ///< Moveconstructor.

    ~OrderError() noexcept override {} ///< Destructor (empty).

private:
    void operator=(const OrderError&) = delete;
    void operator=(OrderError&&) = delete;
};

} // End namespace newton_ns

// ************************************************************************** //

#endif // EXCEPTIONS_H

