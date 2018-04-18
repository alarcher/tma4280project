#ifndef QUADRATURE_QS1K07_H_
#define QUADRATURE_QS1K07_H_

namespace dolfin
{

//-----------------------------------------------------------------------------
typedef quadrature<interval, 7, 4> QS1K07;
//-----------------------------------------------------------------------------
template<> inline real const * QS1K07::weights()
{
  static real const _v[] =
    { + 0.1739274225687269,
      + 0.3260725774312730,
      + 0.3260725774312730,
      + 0.1739274225687269 };
  return _v;
}
//-----------------------------------------------------------------------------
template<> inline real const * QS1K07::points()
{
  static real const _v[] =
    { 0.0694318442029737,
      0.3300094782075719,
      0.6699905217924281,
      0.9305681557970262 };
  return _v;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* QUADRATURE_QS1K07_H_ */
