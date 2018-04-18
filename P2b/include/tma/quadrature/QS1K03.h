#ifndef QUADRATURE_QS1K03_H_
#define QUADRATURE_QS1K03_H_

namespace dolfin
{

//-----------------------------------------------------------------------------
typedef quadrature<interval, 3, 2> QS1K03;
//-----------------------------------------------------------------------------
template<> inline real const * QS1K03::weights()
{
  static real const _v[] =
    { +0.500000000000000,
      +0.500000000000000 };
  return _v;
}
//-----------------------------------------------------------------------------
template<> inline real const * QS1K03::points()
{
  static real const _v[] =
    { 0.2113248654051871,
      0.7886751345948128 };
  return _v;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* QUADRATURE_QS1K03_H_ */
