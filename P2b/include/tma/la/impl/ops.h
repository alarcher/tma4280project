#ifndef TMA_LA_MATRIX_OPS_H_
#define TMA_LA_MATRIX_OPS_H_

namespace tma
{

/******************************************************************************
 * Operators
 ******************************************************************************/

template<class A, class B, A&(*O)(A&, B const&)>
struct __la_unary_op
{
  __la_unary_op(B const& b) : _b(b) { }
  inline A& operator()(A& a) { return O(a, _b); }
  B const&  _b;
};

template<class A, class B, class C, A&(*O)(A&, B const&, C const&)>
struct __la_binary_op
{
  __la_binary_op(B const& b, C const& c) : _b(b), _c(c) { }
  inline A& operator()(A& a) { return O(a, _b, _c); }
  B const&  _b;
  C const&  _c;
};

template<class A, class B, class C, class D, A&(*O)(A&, B const&, C const&, D const&)>
struct __la_ternary_op
{
  __la_ternary_op(B const& b, C const& c, D const& d) : _b(b), _c(c), _d(d) { }
  inline A& operator()(A& a) { return O(a, _b, _c, _d); }
  B const&  _b;
  C const&  _c;
  D const&  _d;
};

//-----------------------------------------------------------------------------

} /* namespace */

#endif /* TMA_LA_MATRIX_H_ */
