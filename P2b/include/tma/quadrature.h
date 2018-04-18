#ifndef TMA_QUADRATURE_H_
#define TMA_QUADRATURE_H_

#include <tma/types.h>

namespace tma
{

template<class T> struct QR;

} /* namespace tma */

#include <tma/cell/interval.h>

namespace tma
{

template<>
struct QR<interval>
{

};

} /* namespace tma */

#include <tma/cell/triangle.h>

namespace tma
{

template<>
struct QR<triangle>
{

};

//-----------------------------------------------------------------------------
template<class Shape, uint K, uint N>
struct quadrature
{
  static uint const d = Shape::topological_dimension;
  static uint const k = K;
  static uint const n = N;
  //---
  static real const * weights() { return NULL; }
  static real const * points()  { return NULL; }
};

//-----------------------------------------------------------------------------
// interval
#include "quadrature/QS1K01.h" // 1 point
#include "quadrature/QS1K03.h" // 2 points
#include "quadrature/QS1K05.h" // 3 points
#include "quadrature/QS1K07.h" // 4 points
#include "quadrature/QS1K09.h" // 5 points
//-----------------------------------------------------------------------------
// triangle
#include "quadrature/QS2K01.h" // 1 ZT 1 point
#include "quadrature/QS2K02.h" // 2 SF 3 points
#include "quadrature/QS2K03.h" // 3 SF 4 points
#include "quadrature/QS2K04.h" // 4 SF 6 points
#include "quadrature/QS2K05.h" // 5 SF 7 points
#include "quadrature/QS2K06.h" // 6 SF 12 points
#include "quadrature/QS2K07.h" // 7    13 points
//-----------------------------------------------------------------------------

} /* namespace tma */

#endif /* TMA_QUADRATURE_H_ */

