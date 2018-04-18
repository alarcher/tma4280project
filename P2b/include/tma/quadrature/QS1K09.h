#ifndef QUADRATURE_QS1K09_H_
#define QUADRATURE_QS1K09_H_

#include <cmath>

namespace dolfin
{

real const a = 1./3.*std::sqrt(5.+4.*std::sqrt(5./14.)) ;
real const b = 1./3.*std::sqrt(5.-4.*std::sqrt(5./14.)) ;
real const wa = (161./450.-13./(180.*std::sqrt(5./14.)))/2.;
real const wb = (161./450.+13./(180.*std::sqrt(5./14.)))/2. ;

//-----------------------------------------------------------------------------
typedef quadrature<interval, 9, 5> QS1K09;
typedef QS1K09 QS1MAX;
//-----------------------------------------------------------------------------
template<> inline real const * QS1K09::weights()
{
  static real const _v[] =
    { +wa,
      +wb,
      +64.0/225.0,
      +wb,
      +wa };
  return _v;
}
//-----------------------------------------------------------------------------
template<> inline real const * QS1K09::points()
{
  static real const _v[] =
    { (1.0 - a)/2.0,
      (1.0 - b)/2.0,
      0.5,
      (1.0 + b)/2.0,
      (1.0 + a)/2.0 };
  return _v;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* QUADRATURE_QS1K09_H_ */
