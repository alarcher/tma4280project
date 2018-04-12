#ifndef TMA_LA_VECTOR_H_
#define TMA_LA_VECTOR_H_

#include <tma/carray.h>

#include <cmath>
#include <algorithm>

namespace tma
{

/******************************************************************************
 * Structures for binary ops
 ******************************************************************************/

//--- Scalar-Vector Product ---------------------------------------------------
template<class K, class E>
struct __svp_op
{
  __svp_op(K const& a, E const& x) :
      a_(a),
      x_(x)
  {
  }
  inline E& operator()(E& y)
  {
    uint const n = y.dim();
#if LA_ALL_OPENMP
#pragma omp parallel for schedule(static) // For course demo
#endif /* LA_ALL_OPENMP */
    for (uint i = 0; i < n; ++i)
    {
      y[i] = a_ * x_[i];
    }
    return y;
  }
  K const&  a_;
  E const&  x_;
};

//--- Vector-Vector Addition --------------------------------------------------
template<class K, class E>
struct __vva_op
{
  __vva_op(E const& x, E const& y) :
      x_(x),
      y_(y)
  {
  }
  inline E& operator()(E& z)
  {
    uint const n = z.dim();
#if LA_ALL_OPENMP
#pragma omp parallel for schedule(static) // For course demo
#endif /* LA_ALL_OPENMP */
    for (uint i = 0; i < n; ++i)
    {
      z[i] = x_[i] + y_[i];
    }
    return z;
  }
  E const& x_;
  E const& y_;
};

/******************************************************************************
 * Vector implementation
 ******************************************************************************/

template<class K>
struct __dense_vector
{
  __dense_vector(uidx n) :
      n_(n),
      V_(zcarray<K>(n))
  {
  }

  ~__dense_vector()
  {
    delete[] V_;
  }

  inline K& operator[](uidx i)
  {
    return V_[i];
  }
  inline K operator[](uidx i) const
  {
    return V_[i];
  }

  inline uidx dim(uint i = 0) const
  {
   return n_;
  }

  inline __dense_vector <K>& operator=(__dense_vector <K> const& other)
  {
    std::copy(other.V_, other.V_ + n_, V_);
    return *this;
  }

  inline __dense_vector <K>& operator=(K a)
  {
    std::fill(V_, V_ + n_, a);
    return *this;
  }

  inline __dense_vector <K>& operator+=(K a)
  {
    for (K * v = V_; v != V_ + n_; ++v) { *v += a; }
    return *this;
  }

  inline __dense_vector <K>& operator-=(K a)
  {
    for (K * v = V_; v != V_ + n_; ++v) { *v -= a; }
    return *this;
  }

  inline __dense_vector <K>& operator*=(K a)
  {
    for (K * v = V_; v != V_ + n_; ++v) { *v *= a; }
    return *this;
  }

  inline __dense_vector <K>& operator/=(K a)
  {
    K const b = K(1)/a;
    for (K * v = V_; v != V_ + n_; ++v) { *v *= b; }
    return *this;
  }

  template<class BinaryOp>
  inline __dense_vector <K>& operator=(BinaryOp op)
  {
    return op(*this);
  }

private:

  uidx const n_;
  K * const V_;
};

/******************************************************************************
 * Operations
 ******************************************************************************/

//--- Inner Product -----------------------------------------------------------

template<class K>
inline K
operator*(__dense_vector <K> const& x, __dense_vector <K> const& y)
{
  uidx const n = y.dim();
  K val = 0.0;
#if LA_ALL_OPENMP
#pragma omp parallel for schedule(static) reduction(+:val) // For course demo
#endif /* LA_ALL_OPENMP */
  for (uint i = 0; i < n; ++i)
  {
    val += x[i] * y[i];
  }
  return val;
}

//--- Scalar-Vector Product ---------------------------------------------------

template<class K>
inline __svp_op <K, __dense_vector <K> >
operator*(K const& a, __dense_vector <K> const& x)
{
  return __svp_op <K, __dense_vector <K> >(a, x);
}

//--- Vector-Vector Addition --------------------------------------------------

template<class K>
inline __vva_op <K, __dense_vector <K> >
operator+(__dense_vector <K> const& x, __dense_vector <K> const& y)
{
  return __vva_op<K, __dense_vector <K> >(x, y);
}

//-----------------------------------------------------------------------------

#define vector  __dense_vector

//-----------------------------------------------------------------------------

} /* namespace */

#endif /* TMA_LA_VECTOR_H_ */
