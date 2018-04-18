#ifndef QUADRATURE_QS1K01_H_
#define QUADRATURE_QS1K01_H_

namespace dolfin
{

//-----------------------------------------------------------------------------
typedef quadrature<interval, 1, 1> QS1K01;
//-----------------------------------------------------------------------------
template<> inline real const * QS1K01::weights()
{
  static real const _v[] =
    { + 1.0 };
  return _v;
}
//-----------------------------------------------------------------------------
template<> inline real const * QS1K01::points()
{
  static real const _v[] =
    { + 0.5 };
  return _v;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* QUADRATURE_QS1K01_H_ */
