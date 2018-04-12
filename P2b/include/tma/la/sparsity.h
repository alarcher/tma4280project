#ifndef TMA_LA_SPARSITY_H_
#define TMA_LA_SPARSITY_H_

#include <tma/la/vector.h>

#include <set>

namespace tma
{

/******************************************************************************
 * Sparsity pattern
 ******************************************************************************/

struct sparsity_pattern
{
  sparsity_pattern(uidx m, uidx n) :
    m_(m),
    n_(n),
    z_(0),
    p_(m ? new std::set<uidx>[m] : NULL)
  {
  }

  ~sparsity_pattern()
  {
    delete [] p_;
  }

  uidx dim(uint i) const { return i ? n_ : m_; }

  void insert(uidx i, uidx const * c0, uidx const * c1)
  {
    uidx n = p_[i].size();
    p_[i].insert(c0, c1);
    z_ += p_[i].size() - n;
  }

  void get(uidx i, uidx * c1) const
  {
    std::copy(p_[i].begin(), p_[i].end(), c1);
  }

  uidx nz(uint i) const { return p_[i].size(); }

  uidx nz() const { return z_; }

private:

  uidx const m_;
  uidx const n_;
  uidx z_;
  std::set<uidx> * p_;
};

} /* namespace */

#endif /* TMA_LA_SPARSITY_H_ */

