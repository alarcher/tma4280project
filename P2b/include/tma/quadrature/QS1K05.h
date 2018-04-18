#ifndef QUADRATURE_QS1K05_H_
#define QUADRATURE_QS1K05_H_

namespace dolfin
{

//-----------------------------------------------------------------------------
typedef quadrature<interval, 5, 3> QS1K05;
//-----------------------------------------------------------------------------
template<> inline real const * QS1K05::weights()
{
  static real const _v[] =
    { + 0.277777777777778,
      + 0.444444444444444,
      + 0.277777777777778 };
  return _v;
}
//-----------------------------------------------------------------------------
template<> inline real const * QS1K05::points()
{
  static real const _v[] =
    { 0.112701665379258,
      0.500000000000000,
      0.887298334620742 };
  return _v;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* QUADRATURE_QS1K05_H_ */
