#ifndef TMA_LA_MATRIX_DENSE_R_H_
#define TMA_LA_MATRIX_DENSE_R_H_

#include <tma/la/vector.h>

namespace tma
{

/******************************************************************************
 * Dense row-major matrix
 ******************************************************************************/

template<class K>
struct __dense_rmatrix
{
  __dense_rmatrix(uidx m, uidx n) :
      m_(m),
      n_(n),
      V_(m * n),
      M_(m * n ? ncarray<K*>(m) : NULL)
  {
    if (M_)
    {
      M_[0] = &V_[0];
      for (uidx i = 0; i < m - 1; ++i) { M_[i + 1] = M_[i] + n; }
    }
  }

  ~__dense_rmatrix()
  {
    delete[] M_;
  }

  // Row vectors
  inline K* operator[](uidx i)
  {
    return M_[i];
  }
  inline K const* operator[](uidx i) const
  {
    return M_[i];
  }

  // Entry access
  inline K& operator()(uidx i, uidx j)
  {
    return M_[i][j];
  }

  inline K operator()(uidx i, uidx j) const
  {
    return M_[i][j];
  }

  inline uidx dim(uint i) const
  {
   return i ? n_ : m_;
  }

  explicit operator vector<K>&()
  {
    return V_;
  }

  explicit operator vector<K> const&() const
  {
    return V_;
  }

  inline __dense_rmatrix <K>& operator=(__dense_rmatrix <K> const& other)
  {
    std::copy(other.V_, other.V_ + m_ * n_, V_);
    return *this;
  }

  inline __dense_rmatrix <K>& operator=(K a)
  {
    std::fill(V_, V_ + m_ * n_, a);
    return *this;
  }

private:

  uidx const m_;
  uidx const n_;
  vector<K>  V_;
  K ** const M_;

  //---------------------------------------------------------------------------
  friend struct __mvp_op<K, __dense_rmatrix <K>, vector<K> >;
  inline void mvp(vector<K> const& x, vector<K>& y) const
  {
    for (uint i = 0; i < m_; ++i)
    {
      y[i] = 0.;
      for (uint j = 0; j < n_; ++j)
      {
        y[i] += M_[i][j] * x[j];
      }
    }
  }

  //---------------------------------------------------------------------------
};

/******************************************************************************
 * Operations
 ******************************************************************************/

//--- Inner Product -----------------------------------------------------------

template<class K>
inline K
operator&(__dense_rmatrix <K> const& A, __dense_rmatrix <K> const& B)
{
  return static_cast<vector <K> const&>(A) * static_cast<vector <K> const&>(B);
}

//--- Matrix-Vector Product ---------------------------------------------------

template<class K>
inline __mvp_op <K, __dense_rmatrix <K>, vector <K> >
operator*(__dense_rmatrix <K> const& A, vector <K> const& x)
{
  return __mvp_op <K,  __dense_rmatrix <K>, vector <K> >(A , x);
}

//-----------------------------------------------------------------------------

} /* namespace */

#endif /* TMA_LA_MATRIX_DENSE_R_H_ */
