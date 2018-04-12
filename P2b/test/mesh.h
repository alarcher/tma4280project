#ifndef TMA_TEST_MESH_H_
#define TMA_TEST_MESH_H_

#include <tma.h>

namespace tma
{

namespace tests
{

struct mesh
{
  static void test()
  {
    // interval mesh
    message("mesh::interval");
    {
      // Interval with two subintervals
      topology<interval> T(2, 3);
      T(0)[0] = 0;
      T(0)[1] = 1;
      T(1)[0] = 1;
      T(1)[1] = 2;
      T.dump();
      geometry<1> G(3);
      G(0)[0] = 0.0;
      G(1)[0] = 0.5;
      G(2)[0] = 1.0;
      G.dump();

      // Test mesh constructor
      tma::mesh<interval, 1> (2, 3);

      // Unit interval
      unit_interval M(4);
      M.dump();

      // Affine mapping
      AffineMapping<interval, 1> F;
      for (uint i = 0; i < M.ncells(); ++i)
      {
        std::cout << i<< ":\n";
        F(M, i); F.disp();
      }
    }

    // triangle mesh
    message("mesh::triangle");
    {
      // A square with two triangles would use this topology
      topology<triangle> T(2, 4);
      T(0)[0] = 0;
      T(0)[1] = 1;
      T(0)[2] = 2;
      T(1)[0] = 0;
      T(1)[1] = 2;
      T(1)[2] = 3;
      T.dump();

      geometry<2> G(4);
      G(0)[0] = 0.0; G(0)[1] = 0.0;
      G(1)[0] = 1.0; G(1)[1] = 0.0;
      G(2)[0] = 1.0; G(2)[1] = 1.0;
      G(3)[0] = 0.0; G(3)[1] = 1.0;
      G.dump();

      // Test mesh constructor
      tma::mesh<triangle, 2> (2, 4);

      // Unit interval
      unit_square M(4, 2);
      M.dump();

      // Affine mapping
      AffineMapping<triangle, 2> F;
      for (uint i = 0; i < M.ncells(); ++i)
      {
        std::cout << i<< ":\n";
        F(M, i); F.disp();
      }
    }
  }
};

} /* namespace tests */

} /* namespace tma */

#endif /* TMA_TEST_MESH_H_ */
