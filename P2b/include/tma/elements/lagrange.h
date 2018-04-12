#ifndef TMA_ELEMENTS_LAGRANGE_H_
#define TMA_ELEMENTS_LAGRANGE_H_

#include <tma/types.h>

#include <tma/mesh/mapping.h>

namespace tma
{

template<class T> struct P1;

} /* namespace tma */

#include <tma/cell/interval.h>

namespace tma
{

template<>
struct P1<interval>
{
  /// Return the dimension of the reference element
  uint dim() const
  {
    return 2;
  }

  // Coordinates of degrees of freedom in the reference interval
  real const * x(uint i)
  {
    static real const _s[2][1] = { { 0.0 }, { 1.0 } };
    return _s[i];
  }

  // Evaluate the polynomial basis in the reference element
  void operator()(real const* x, real* v) const
  {
    v[0] = 1.0 - x[0];
    v[1] = x[0];
  }

  // Evaluate the polynomial basis derivatives in the reference element
  void d(real const* x, real(& v)[2][1]) const
  {
    v[0][0] = - 1.0;
    v[1][0] = + 1.0;
  }

  // Evaluate the polynomial basis derivatives in the reference element (flat)
  void d(real const* x, real* v) const
  {
    v[0] = - 1.0;
    v[1] = + 1.0;
  }

  // Local equation
  struct local
  {
    uidx basis[2];
    uidx shape[2];
    real matrix[2][2];

    void operator()(mesh<interval, 1> const& m, uint i)
    {
      std::copy(m.K(i), m.K(i) + 2, basis);
      std::copy(m.K(i), m.K(i) + 2, shape);
    }
  };

  typedef AffineMapping<interval, 1> mapping;
};

} /* namespace tma */

#include <tma/cell/triangle.h>

namespace tma
{

template<>
struct P1<triangle>
{
  /// Return the dimension of the reference element
  uint dim() const
  {
    return 3;
  }

  // Coordinates of degrees of freedom in the reference triangle
  real const * x(uint i)
  {
    static real const x_[3][2] = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 } };
    return x_[i];
  }

  // Evaluate the polynomial basis at x
  void operator()(real const* x, real* v)
  {
    v[0] = 1.0 - x[0] - x[1];
    v[1] = x[0];
    v[2] = x[1];
  }

  // Evaluate the polynomial basis derivatives in the reference element
  void d(real const* x, real(& v)[3][2])
  {
    v[0][0] = -1.0;
    v[0][1] = -1.0;
    v[1][0] = +1.0;
    v[1][1] =  0.0;
    v[2][0] =  0.0;
    v[2][1] = +1.0;
  }

  // Evaluate the polynomial basis derivatives in the reference element (flat)
  void d(real const* x, real* v)
  {
    v[0] = -1.0;
    v[1] = -1.0;
    v[2] = +1.0;
    v[3] =  0.0;
    v[4] =  0.0;
    v[5] = +1.0;
  }

  // Local equation
  struct local
  {
    uidx basis[3];
    uidx shape[3];
    real matrix[3][3];

    void operator()(mesh<triangle, 2> const& m, uint i)
    {
      std::copy(m.K(i), m.K(i) + 3, basis);
      std::copy(m.K(i), m.K(i) + 3, shape);
    }
  };

  typedef AffineMapping<triangle, 2> mapping;
};

} /* namespace tma */

#endif /* TMA_ELEMENTS_LAGRANGE_H_ */
