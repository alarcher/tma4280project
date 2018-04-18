#ifndef QUADRATURE_QS3K01_H_
#define QUADRATURE_QS3K01_H_

namespace dolfin
{

//-----------------------------------------------------------------------------
typedef quadrature<tetrahedron, 1, 1> QS3K01;
//-----------------------------------------------------------------------------
template<> inline real const * QS3K01::weights()
{
  static real const _v[] =
    { + 0.166666666667 };
  return _v;
}
//-----------------------------------------------------------------------------
template<> inline real const * QS3K01::points()
{
  static real const _v[] =
    { 0.250000000000, 0.250000000000, 0.250000000000 };
  return _v;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* QUADRATURE_QS3K01_H_ */
