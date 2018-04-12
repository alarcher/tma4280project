#ifndef TMA_TEST_LA_H_
#define TMA_TEST_LA_H_

#include <tma.h>

namespace tma
{

namespace tests
{

template<class T>
void test_matrix(T& A)
{
  A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
  A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;
  A(2, 0) = 7; A(2, 1) = 8; A(2, 2) = 9;
  message("matrix created");
  message("%e %e %e", A(0, 0), A(0, 1), A(0, 2));
  message("%e %e %e", A(1, 0), A(1, 1), A(1, 2));
  message("%e %e %e", A(2, 0), A(2, 1), A(2, 2));
  vector<real>  x(3);
  x = 1.0;
  vector<real>  y(3);
  y = A * x;
  message("y =  A x  = [ %e, %e, %e ]^T", y[0], y[1], y[2]);
  tma_assert(y[0] ==  6.);
  tma_assert(y[1] == 15.);
  tma_assert(y[2] == 24.);
}

struct la
{
  static void test()
  {
    // real vector
    {
      vector<real> x(3);
      x[0] = 0.;
      x[1] = 1.;
      x[2] = 2.;
      vector<real> y(3);
      y[0] = 0.;
      y[1] = 2.;
      y[2] = 4.;
      real alpha = x * y;
      // Binary rationals should be exact
      message("alpha = x . y = %e", alpha);
      tma_assert(alpha == 10.);
      vector<real> z(3);
      z = x + y;
      message("z = x + y = [ %e, %e, %e ]^T", z[0], z[1], z[2]);
      tma_assert(z[0] == 0.);
      tma_assert(z[1] == 3.);
      tma_assert(z[2] == 6.);
      z *= 2.0;
      message("z = 2 * z = [ %e, %e, %e ]^T", z[0], z[1], z[2]);
      tma_assert(z[0] ==  0.);
      tma_assert(z[1] ==  6.);
      tma_assert(z[2] == 12.);
      z /= 2.0;
      message("z = z / 2 = [ %e, %e, %e ]^T", z[0], z[1], z[2]);
      tma_assert(z[0] == 0.);
      tma_assert(z[1] == 3.);
      tma_assert(z[2] == 6.);
      z += 2.0;
      message("z = z + 2 = [ %e, %e, %e ]^T", z[0], z[1], z[2]);
      tma_assert(z[0] == 2.);
      tma_assert(z[1] == 5.);
      tma_assert(z[2] == 8.);
      z -= 2.0;
      message("z = z + 2 = [ %e, %e, %e ]^T", z[0], z[1], z[2]);
      tma_assert(z[0] == 0.);
      tma_assert(z[1] == 3.);
      tma_assert(z[2] == 6.);
    }
    // dense real matrix, column-major
    {
      cmatrix<real> A(3, 3);
      test_matrix(A);
    }
    // dense real matrix, row-major
    {
      rmatrix<real> A(3, 3);
      test_matrix(A);
    }
    // sparse real matrix, row-major (used as dense)
    {
      uidx c[9] = { 0, 1, 2, 0, 1, 2, 0, 1, 2 };
      smatrix<real> A(3, 3, 3, c);
      tma_assert(A.nz()    ==  9);
      tma_assert(A.nzmin() ==  3);
      tma_assert(A.nzmax() ==  3);
      for (uint i = 0; i < A.dim(0); ++i) { tma_assert(A.nz(i) ==  3); }
      test_matrix(A);
    }
    // sparse real matrix, row-major (used as dense)
    {
      uidx z[13] = { 9, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2 };
      smatrix<real> A(3, 3, z);
      tma_assert(A.nz()    ==  9);
      tma_assert(A.nzmin() ==  3);
      tma_assert(A.nzmax() ==  3);
      for (uint i = 0; i < A.dim(0); ++i) { tma_assert(A.nz(i) ==  3); }
      test_matrix(A);
    }
    // sparse real matrix, row-major (identity)
    {
      uidx z[7] = { 3, 1, 0, 1, 1, 1, 2 };
      smatrix<real> A(3, 3, z);
      tma_assert(A.nz()    ==  3);
      tma_assert(A.nzmin() ==  1);
      tma_assert(A.nzmax() ==  1);
      for (uint i = 0; i < A.dim(0); ++i) { tma_assert(A.nz(i) ==  1); }
      A(0,0) = 1.0; A(1,1) = 1.0; A(2,2) = 1.0;
      rmatrix<real> B(3, 3);
      B(0, 0) = 1; B(0, 1) = 2; B(0, 2) = 3;
      B(1, 0) = 4; B(1, 1) = 5; B(1, 2) = 6;
      B(2, 0) = 7; B(2, 1) = 8; B(2, 2) = 9;
      rmatrix<real> C(3, 3);
      //C = A * B;
    }
    // create sparsity pattern from vertex connectivities in 1d
    {
      unit_interval M(4);
      sparsity_pattern sp(M.nverts(), M.nverts());
      for (uint i = 0; i < M.ncells(); ++i)
      {
        for (uint v = 0; v < interval::num(0); ++v)
        {
          sp.insert(M.K(i)[v], M.K(i), M.K(i) + interval::num(0));
        }
      }
    }
    // create sparsity pattern from vertex connectivities in 2d
    {
      unit_square M(4, 4);
      sparsity_pattern sp(M.nverts(), M.nverts());
      for (uint i = 0; i < M.ncells(); ++i)
      {
        for (uint v = 0; v < triangle::num(0); ++v)
        {
          sp.insert(M.K(i)[v], M.K(i), M.K(i) + triangle::num(0));
        }
      }
    }
  }
};

} /* namespace tests */

} /* namespace tma */

#endif /* TMA_TEST_LA_H_ */
