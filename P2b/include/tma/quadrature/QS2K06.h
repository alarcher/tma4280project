#ifndef QUADRATURE_QS2K06_H_
#define QUADRATURE_QS2K06_H_

namespace dolfin
{

//-----------------------------------------------------------------------------
typedef quadrature<triangle, 6, 12> QS2K06;
//-----------------------------------------------------------------------------
template<> inline real const * QS2K06::weights()
{
  static real const _v[] =
    { + 0.0583931378631895,
      + 0.0583931378631895,
      + 0.0583931378631895,
      + 0.0254224531851035,
      + 0.0254224531851035,
      + 0.0254224531851035,
      + 0.0414255378091870,
      + 0.0414255378091870,
      + 0.0414255378091870,
      + 0.0414255378091870,
      + 0.0414255378091870,
      + 0.0414255378091870 };
  return _v;
}
//-----------------------------------------------------------------------------
template<> inline real const * QS2K06::points()
{
  static real const _v[] =
    { 0.501426509658179, 0.249286745170910,
      0.249286745170910, 0.249286745170910,
      0.249286745170910, 0.501426509658179,
      0.873821971016996, 0.063089014491502,
      0.063089014491502, 0.063089014491502,
      0.063089014491502, 0.873821971016996,
      0.053145049844817, 0.310352451033784,
      0.310352451033784, 0.053145049844817,
      0.053145049844817, 0.636502499121399,
      0.636502499121399, 0.053145049844817,
      0.636502499121399, 0.310352451033784,
      0.310352451033784, 0.636502499121399 };
  return _v;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* QUADRATURE_QS2K06_H_ */
