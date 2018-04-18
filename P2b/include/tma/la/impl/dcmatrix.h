#ifndef TMA_LA_MATRIX_DENSE_C_H_
#define TMA_LA_MATRIX_DENSE_C_H_

#include <tma/la/vector.h>

namespace tma
{

/******************************************************************************
 * Dense column-major matrix
 ******************************************************************************/

template<class K>
struct __dense_cmatrix
{
  __dense_cmatrix(uidx m, uidx n) :
      m_(m),
      n_(n),
      V_(m * n),
      M_(m * n ? ncarray<K*>(n) : NULL)
  {
    if (M_)
    {
      M_[0] = &V_[0];
      for (uidx i = 0; i < n - 1; ++i) { M_[i + 1] = M_[i] + m; }
    }
  }

  ~__dense_cmatrix()
  {
    delete[] M_;
  }

  // Column vectors
  inline K* operator[](uidx j)
  {
    return M_[j];
  }
  inline K const* operator[](uidx j) const
  {
    return M_[j];
  }

  // Entry access
  inline K& operator()(uidx i, uidx j)
  {
    return M_[j][i];
  }

  inline K operator()(uidx i, uidx j) const
  {
    return M_[j][i];
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

  inline __dense_cmatrix <K>& operator=(__dense_cmatrix <K> const& other)
  {
    std::copy(other.V_, other.V_ + m_ * n_, V_);
    return *this;
  }

  inline __dense_cmatrix <K>& operator=(K a)
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
  friend struct __mvp_op<K, __dense_cmatrix <K>, vector<K> >;
  inline void mvp(vector<K> const& x, vector<K>& y) const
  {
    y = 0.;
    for (uint j = 0; j < n_; ++j)
    {
      for (uint i = 0; i < m_; ++i)
      {
        y[i] += M_[j][i] * x[j];
      }
    }
  }

  //---------------------------------------------------------------------------
  friend struct __mmp_op<K, __dense_cmatrix <K> >;
  inline void mmp(__dense_cmatrix <K> const& Q, __dense_cmatrix <K>& A) const
  {
    A = 0.;
    for (uint j = 0; j < Q.n_; ++j)
    {
      for (uint k = 0; k < n_; ++k)
      {
        for (uint i = 0; i < m_; ++i)
        {
          A.M_[j][i] += Q.M_[j][k] * M_[k][i];
        }
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
operator&(__dense_cmatrix <K> const& A, __dense_cmatrix <K> const& B)
{
  return static_cast<vector <K> const&>(A) * static_cast<vector <K> const&>(B);
}

//--- Matrix-Vector Product ---------------------------------------------------

template<class K>
inline __mvp_op <K, __dense_cmatrix <K>, vector <K> >
operator*(__dense_cmatrix <K> const& A, vector <K> const& x)
{
  return __mvp_op <K, __dense_cmatrix <K>, vector <K> >(A , x);
}

//-----------------------------------------------------------------------------

} /* namespace */

#endif /* TMA_LA_MATRIX_DENSE_C_H_ */
