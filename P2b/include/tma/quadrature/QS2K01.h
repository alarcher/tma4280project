#ifndef QUADRATURE_QS2K01_H_
#define QUADRATURE_QS2K01_H_

namespace dolfin
{

//-----------------------------------------------------------------------------
typedef quadrature<triangle, 1, 1> QS2K01;
//-----------------------------------------------------------------------------
template<> inline real const * QS2K01::weights()
{
  static real const _v[] =
    { + 0.5000000000000000000 };
  return _v;
}
//-----------------------------------------------------------------------------
template<> inline real const * QS2K01::points()
{
  static real const _v[] =
    { 0.3333333333333333333, 0.3333333333333333333 };
  return _v;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* QUADRATURE_QS2K01_H_ */
