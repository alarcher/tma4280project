#ifndef QUADRATURE_QS3K03_H_
#define QUADRATURE_QS3K03_H_

namespace dolfin
{

//-----------------------------------------------------------------------------
typedef quadrature<tetrahedron, 3, 5> QS3K03;
//-----------------------------------------------------------------------------
template<> inline real const * QS3K03::weights()
{
  static real const _v[] =
    { - 0.13333333333333333333,
      + 0.07500000000000000000,
      + 0.07500000000000000000,
      + 0.07500000000000000000,
      + 0.07500000000000000000 };
  return _v;
}
//-----------------------------------------------------------------------------
template<> inline real const * QS3K03::points()
{
  static real const _v[] =
    { 0.2500000000000000, 0.2500000000000000, 0.2500000000000000,
      0.1666666666666667, 0.1666666666666667, 0.1666666666666667,
      0.5000000000000000, 0.1666666666666667, 0.1666666666666667,
      0.1666666666666667, 0.5000000000000000, 0.1666666666666667,
      0.1666666666666667, 0.1666666666666667, 0.5000000000000000 };
  return _v;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* QUADRATURE_QS3K03_H_ */
