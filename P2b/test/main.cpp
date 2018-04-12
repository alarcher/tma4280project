
#include <tma.h>

#include "cell.h"
#include "mesh.h"
#include "elements.h"
#include "la.h"

using namespace tma;

int main(int argc, char**argv)
{
  //
  tma::tests::cell::test();

  //
  tma::tests::mesh::test();

  //
  tma::tests::elements::test();

  //
  tma::tests::la::test();

  return 0;
}
