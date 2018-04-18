#ifndef TMA_LA_MATRIX_SPARSE_R_H_
#define TMA_LA_MATRIX_SPARSE_R_H_

#include <tma/la/vector.h>
#include <tma/la/sparsity.h>

#include <iomanip>

namespace tma
{

/******************************************************************************
 * Sparse row-major matrix
 ******************************************************************************/

template<class K>
struct __sparse_rmatrix
{
  // With constant non-zero entries per row
  __sparse_rmatrix(uidx m, uidx n, uint z, uidx c[]) :
      m_(m),
      n_(n),
      z_(m * z),
      zmin_(z),
      zmax_(z),
      V_(z_),
      M_(z_ ? ncarray<K*>(m + 1) : NULL), // Use last entry as bound
      I_(z_ ? ncarray<uidx>(z_)  : NULL)
  {
    if (M_)
    {
      M_[0] = &V_[0];
      for (uidx i = 0; i < m; ++i) { M_[i + 1] = M_[i] + z; }
      std::copy(c, c + z_, I_);
    }
  }

  // With given sparsity pattern: z[0] = total nz + m x ( z, {c} )
  __sparse_rmatrix(uidx m, uidx n, uidx z[]) :
      m_(m),
      n_(n),
      z_(*z++),
      zmin_(0),
      zmax_(0),
      V_(z_),
      M_(z_ ? ncarray<K*>(m + 1) : NULL), // Use last entry as bound
      I_(z_ ? ncarray<uidx>(z_)  : NULL)
  {
    if (M_)
    {
      M_[0] = &V_[0];
      uidx * ci = I_;
      zmin_ = z_; zmax_ = 0;
      for (uidx i = 0; i < m; ++i)
      {
        uint nc = *z++;
        zmin_ = std::min(zmin_, nc);
        zmax_ = std::max(zmax_, nc);
        for (uint j = 0; j < nc; ++j)  *ci++ = *z++;
        M_[i + 1] = M_[i] + nc;
      }
    }
  }

  // With given sparsity pattern
  __sparse_rmatrix(sparsity_pattern const& sp) :
      m_(sp.dim(0)),
      n_(sp.dim(1)),
      z_(sp.nz()),
      zmin_(0),
      zmax_(0),
      V_(z_),
      M_(z_ ? ncarray<K*>(m_ + 1) : NULL), // Use last entry as bound
      I_(z_ ? ncarray<uidx>(z_)  : NULL)
  {
    if (M_)
    {
      M_[0] = &V_[0];
      uidx * ci = I_;
      zmin_ = z_; zmax_ = 0;
      for (uidx i = 0; i < m_; ++i)
      {
        uint nc = sp.nz(i);
        zmin_ = std::min(zmin_, nc);
        zmax_ = std::max(zmax_, nc);
        sp.get(i, ci);
        ci += nc;
        M_[i + 1] = M_[i] + nc;
      }
    }
  }

  ~__sparse_rmatrix()
  {
    delete[] M_;
    delete[] I_;
  }

  // Row vectors (non-zero only)
  inline K* operator[](uidx i)
  {
    return M_[i];
  }
  inline K const* operator[](uidx i) const
  {
    return M_[i];
  }

  // Entry access (slow)
  inline K& operator()(uidx i, uidx j)
  {
    uidx * ci = I_ + (M_[i] - M_[0]);
    for(K * x = M_[i]; x != M_[i + 1]; ++x)
    {
      if (*ci++ == j) return *x;
    }
    perror("Invalid column index");
    return M_[i][0];
  }

  inline K operator()(uidx i, uidx j) const
  {
    uidx * ci = I_ + (M_[i] - M_[0]);
    for(K * x = M_[i]; x != M_[i + 1]; ++x)
    {
      if (*ci++ == j) return *x;
    }
    return K(0);
  }

  inline uidx dim(uint i) const
  {
   return i ? n_ : m_;
  }

  inline uidx nz(uidx i) const
  {
   return M_[i + 1] - M_[i];
  }

  inline uidx nzmin() const
  {
   return zmin_;
  }

  inline uidx nzmax() const
  {
   return zmax_;
  }

  inline uidx nz() const
  {
   return z_;
  }

  explicit operator vector<K>&()
  {
    return V_;
  }

  explicit operator vector<K> const&() const
  {
    return V_;
  }

  inline __sparse_rmatrix <K>& operator=(K a)
  {
    std::fill(V_, V_ + z_, a);
    return *this;
  }

  void disp() const
  {
    uidx * ci = I_;
    for (uint i = 0; i < m_; ++i)
    {
      std::cout << std::setw(8) << i << ": ";
      for (K * aij = M_[i]; aij != M_[i + 1]; ++aij)
      {
        std::cout << std::setw(8) << *ci++;
      }
      std::cout << "\n";
    }
  }

  friend
  K operator& (__sparse_rmatrix <K> const& P, __sparse_rmatrix <K> const& Q)
  {
    K a = 0;
    uidx * cp = P.I_;
    uidx * cq = Q.I_;
    for (uint i = 0; i < P.m_; ++i)
    {
      K * qij = Q.M_[i];
      for (K * pij = P.M_[i]; pij != P.M_[i + 1]; ++pij, ++cp)
      {
        for (; qij != Q.M_[i + 1]; ++qij, ++cq)
        {
          if (*cp == *cq) a+= *pij * qij;
        }
        if (qij == Q.M_[i + 1]) break;
      }
    }
    return a;
  }

private:

  uidx const m_;
  uidx const n_;
  uidx const z_;
  uint zmin_;
  uint zmax_;
  vector<K>  V_;
  K ** const M_;
  uidx *     I_;

  //---------------------------------------------------------------------------
  friend struct __mvp_op<K, __sparse_rmatrix <K>, vector<K> >;
  inline void mvp(vector<K> const& x, vector<K>& y) const
  {
    uidx * ci = I_;
    for (uint i = 0; i < m_; ++i)
    {
      y[i] = 0.;
      for (K * aij = M_[i]; aij != M_[i + 1]; ++aij)
      {
        y[i] += *aij * x[*ci++];
      }
    }
  }

  //---------------------------------------------------------------------------
  friend struct __mmp_op<K, __sparse_rmatrix <K> >;
  inline void mmp(__sparse_rmatrix <K> const& Q, __sparse_rmatrix <K>& A) const
  {
    for (uint i = 0; i < m_; ++i)
    {
      for (uint j = 0; j < Q.n_; ++j)
      {
        A.M_[i][j] = 0.;
        for (uint k = 0; k < n_; ++k)
        {
          A.M_[i][j] += M_[i][k] * Q.M_[k][j];
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
operator&(__sparse_rmatrix <K> const& A, __sparse_rmatrix <K> const& B)
{
  return static_cast<vector <K> const&>(A) * static_cast<vector <K> const&>(B);
}

//--- Matrix-Vector Product ---------------------------------------------------

template<class K>
inline __mvp_op <K, __sparse_rmatrix <K>, vector <K> >
operator*(__sparse_rmatrix <K> const& A, vector <K> const& x)
{
  return __mvp_op <K,  __sparse_rmatrix <K>, vector <K> >(A , x);
}

//-----------------------------------------------------------------------------

} /* namespace */

#endif /* TMA_LA_MATRIX_SPARSE_R_H_ */
