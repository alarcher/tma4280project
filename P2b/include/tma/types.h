#ifndef TMA_TYPES_H_
#define TMA_TYPES_H_

#include <algorithm>
#include <cmath>
#include <climits>
#include <cfloat>
#include <stdint.h>

namespace tma
{

typedef unsigned int uint;
typedef size_t       uidx;
typedef double       real;

/// Return strong relative real comparison:
/// ( |(x - y) / x| < eps ) && ( |(x - y) / y| < eps )
static inline bool srelcmp(real x, real y, real eps)
{
  real const d = std::fabs(x - y);
  real const m = std::min(std::fabs(x), std::fabs(y));
  // Take very small positive values of d and m into account
  if (d < std::numeric_limits<real>::min() ||
      m < std::numeric_limits<real>::min())
  {
    return true;
  }
  // Trying to avoid underflow issues in most common cases
  return (
      (m > std::numeric_limits<real>::epsilon()) ?
          (d / m < std::fabs(eps)) :
          (d / std::sqrt(m) < std::fabs(eps) * std::sqrt(m)));
}

} /* namespace tma */

#endif /* TMA_TYPES_H _ */
