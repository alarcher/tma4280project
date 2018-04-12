#ifndef TMA_LA_MATRIX_H_
#define TMA_LA_MATRIX_H_

#include <tma/types.h>

namespace tma
{

/******************************************************************************
 * Structure for binary ops
 ******************************************************************************/

template<class K, class M, class E>
struct __mvp_op
{
  __mvp_op(M const& A, E const& x) : A_(A), x_(x) { }
  inline E& operator()(E& y) { A_.mvp(x_, y); return y; }
  M const&  A_;
  E const&  x_;
};

} /* namespace */

//-----------------------------------------------------------------------------

#include "impl/dcmatrix.h"

#define cmatrix __dense_cmatrix

//-----------------------------------------------------------------------------

#include "impl/drmatrix.h"

#define rmatrix __dense_rmatrix

//-----------------------------------------------------------------------------

#include "impl/srmatrix.h"

#define smatrix __sparse_rmatrix

//-----------------------------------------------------------------------------

#endif /* TMA_LA_MATRIX_H_ */
