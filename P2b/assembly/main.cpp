
#include <tma.h>

using namespace tma;

template<class T, uint D>
void assembly(mesh<T, D>& m)
{
  real t;
  P1<T>                   FE; // Finite element
  typename P1<T>::local   EQ; // Local equation
  typename P1<T>::mapping TM; // Transport mapping

  sparsity_pattern SP(m.nverts(), m.nverts());
  tic(t);
  for (uint c = 0; c < m.ncells(); ++c)
  {
    // Set local equation and compute mapping
    EQ(m, c);
    for (uint i = 0; i < FE.dim(); ++i)
    {
      SP.insert(EQ.basis[i], EQ.shape, EQ.shape + FE.dim());
    }
  }
  tic(t);
  message("Sparsity: %e seconds", t);
  smatrix<real> A(SP);
  A.disp();
  tic(t);
  message("Matrix A: %e seconds", t);
  vector<real>  b(SP.dim(0));
  tic(t);
  message("Vector b: %e seconds", t);
  vector<real>  x(SP.dim(1));
  tic(t);
  message("Vector x: %e seconds", t);
  tic(t);
  for (uint c = 0; c < m.ncells(); ++c)
  {
    // Set local equation
    EQ(m, c);

    // Compute mapping
    TM(m, c);

    // Loop on quadrature points

    // Add to linear system
  }
  tic(t);
  message("Assembly: %e seconds", t);
}

int main(int argc, char**argv)
{
  // interval cell
  message("unit_interval");
  {
    unit_interval m(1 << 20);
    assembly(m);
  }

  // triangle cell
  message("unit_square");
  {
    unit_square m(1 << 10, 1 << 10);
    assembly(m);
  }

  return 0;
}
