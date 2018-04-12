#ifndef TMA_MESH_UNIT_H_
#define TMA_MESH_UNIT_H_

#include <tma/mesh/mesh.h>
#include <tma/cell/interval.h>

namespace tma
{

struct unit_interval : public mesh<interval, 1>
{
  unit_interval(uint n) :
    mesh<interval, 1>(n, n + 1)
  {
    // Topology
    for (uint i = 0; i < n; ++i)
    {
      K(i)[0] = i;
      K(i)[1] = i + 1;
    }
    // Geometry
    for (uint i = 0; i <= n; ++i)
    {
      x(i)[0] = static_cast<real>(i) / n;
    }
  }
};

struct unit_square : public mesh<triangle, 2>
{
  unit_square(uint m, uint n) :
    mesh<triangle, 2>(2*m*n, (m + 1)*(n + 1))
  {
    // Topology
    uint c = 0;
    for (uint j = 0; j < n; ++j)
    {    
      for (uint i = 0; i < m; ++i, ++c)
      {
        uint s = j * (m + 1) + i;
        K(c)[0] = s;
        K(c)[1] = s + 1;
        K(c)[2] = s + m + 2;
        ++c;
        K(c)[0] = s;
        K(c)[1] = s + m + 1;
        K(c)[2] = s + m + 2;
      }
    }
    // Geometry
    uint v = 0;
    for (uint j = 0; j <= n; ++j)
    {
      real y = static_cast<real>(j) / n;
      for (uint i = 0; i <= m; ++i, ++v)
      {
        x(v)[0] = static_cast<real>(i) / m;
        x(v)[1] = y;
      }
    }
  }
};

} /* namespace */


#endif /* TMA_MESH_H_ */
