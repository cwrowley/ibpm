/// \file typedefs.h
/// \brief Header containing template aliases for this package.
/// \author Kevin K. Chen (University of Southern California)

#ifndef TEMPLATE_ALIASES_H
#define TEMPLATE_ALIASES_H

#include <functional>

// ************************************************************************** //

namespace alias {

using std::function;

// ************************************************************************** //

template <typename Vec>
using Function = function<Vec (const Vec&)>;

template <typename Vec>
using InnerProduct = function<double (const Vec&, const Vec&)>;

template <typename Vec>
using Norm = function<double (const Vec&)>;

// ************************************************************************** //

} // End namespace alias

#endif // TEMPLATE_ALIASES_H
