#ifndef TMA_TEST_CELL_H_
#define TMA_TEST_CELL_H_

#include <tma.h>

namespace tma
{

namespace tests
{

struct cell
{
  static void test()
  {
    // interval cell
    message("cell::interval");
    {
      interval C;
      tma_assert(C.dim()  == 1);
      tma_assert(C.num(0) == 2);
      tma_assert(C.num(1) == 1);

      interval::reference R;
      tma_assert(R.x(0)[0] == 0.0);
      tma_assert(R.x(1)[0] == 1.0);
    }

    // triangle cell
    message("cell::triangle");
    {
      triangle C;
      tma_assert(C.dim()  == 2);
      tma_assert(C.num(0) == 3);
      tma_assert(C.num(1) == 3);
      tma_assert(C.num(2) == 1);

      triangle::reference R;
      tma_assert(R.x(0)[0] == 0.0);
      tma_assert(R.x(0)[1] == 0.0);
      tma_assert(R.x(1)[0] == 1.0);
      tma_assert(R.x(1)[1] == 0.0);
      tma_assert(R.x(2)[0] == 0.0);
      tma_assert(R.x(2)[1] == 1.0);
    }
  }
};

} /* namespace tests */

} /* namespace tma */

#endif /* TMA_TEST_CELL_H_ */
