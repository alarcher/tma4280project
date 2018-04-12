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

} /* namespace tma */

#endif /* TMA_QUADRATURE_H_ */

