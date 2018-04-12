#ifndef TMA_TEST_ELEMENTS_H_
#define TMA_TEST_ELEMENTS_H_

#include <tma.h>

namespace tma
{

namespace tests
{

struct elements
{
  static void test()
  {
    // P1 on interval cell
    message("elements::P1<interval>");
    {
      P1<interval> E;
      // Check element dimension
      tma_assert(E.dim()   == 2);
      // Check element node coordinates
      tma_assert(E.x(0)[0] == 0.0);
      tma_assert(E.x(1)[0] == 1.0);
      // Check basis functions
      {
        real v[2];
        E(E.x(0), v);
        tma_assert(v[0]  == 1.0);
        tma_assert(v[1]  == 0.0);
        E(E.x(1), v);
        tma_assert(v[0]  == 0.0);
        tma_assert(v[1]  == 1.0);

        // Check sum
        real x[1];
        for (uint i = 0; i <= 100; ++i)
        {
          x[0] = real(i)/100;
          E(x, v);
          tma_assert(srelcmp(v[0] + v[1], 1.0, 2e-16));
        }
      }
      // Check basis functions derivatives
      {
        real v[2][1];
        E.d(E.x(0), v);
        tma_assert(v[0][0]  == -1.0);
        tma_assert(v[1][0]  == +1.0);
        E.d(E.x(1), v);
        tma_assert(v[0][0]  == -1.0);
        tma_assert(v[1][0]  == +1.0);
      }
    }

    // P1 on triangle cell
    message("elements::P1<triangle>");
    {
      P1<triangle> E;
      // Check element dimension
      tma_assert(E.dim()   == 3);
      // Check element node coordinates
      tma_assert(E.x(0)[0] == 0.0);
      tma_assert(E.x(0)[1] == 0.0);
      tma_assert(E.x(1)[0] == 1.0);
      tma_assert(E.x(1)[1] == 0.0);
      tma_assert(E.x(2)[0] == 0.0);
      tma_assert(E.x(2)[1] == 1.0);
      // Check basis functions
      {
        real v[3];
        E(E.x(0), v);
        tma_assert(v[0]  == 1.0);
        tma_assert(v[1]  == 0.0);
        tma_assert(v[2]  == 0.0);
        E(E.x(1), v);
        tma_assert(v[0]  == 0.0);
        tma_assert(v[1]  == 1.0);
        tma_assert(v[2]  == 0.0);
        E(E.x(2), v);
        tma_assert(v[0]  == 0.0);
        tma_assert(v[1]  == 0.0);
        tma_assert(v[2]  == 1.0);
        // Check sum
        real x[2];
        for (uint i = 0; i <= 100; ++i)
        {
          x[0] = real(i)/100;
          for (uint j = 0; j <= 100 - i; ++j)
          {
            x[1] = real(j)/100;
            E(x, v);
            tma_assert(srelcmp(v[0] + v[1] + v[2], 1.0, 4e-16));
          }
        }
      }
      // Check basis functions derivatives
      {
        real v[3][2];
        E.d(E.x(0), v);
        tma_assert(v[0][0]  == -1.0);
        tma_assert(v[0][1]  == -1.0);
        tma_assert(v[1][0]  == +1.0);
        tma_assert(v[1][1]  ==  0.0);
        tma_assert(v[2][0]  ==  0.0);
        tma_assert(v[2][1]  == +1.0);
        E.d(E.x(1), v);
        tma_assert(v[0][0]  == -1.0);
        tma_assert(v[0][1]  == -1.0);
        tma_assert(v[1][0]  == +1.0);
        tma_assert(v[1][1]  ==  0.0);
        tma_assert(v[2][0]  ==  0.0);
        tma_assert(v[2][1]  == +1.0);
        E.d(E.x(2), v);
        tma_assert(v[0][0]  == -1.0);
        tma_assert(v[0][1]  == -1.0);
        tma_assert(v[1][0]  == +1.0);
        tma_assert(v[1][1]  ==  0.0);
        tma_assert(v[2][0]  ==  0.0);
        tma_assert(v[2][1]  == +1.0);
      }
    }
  }
};

} /* namespace tests */

} /* namespace tma */

#endif /* TMA_TEST_ELEMENTS_H_ */
