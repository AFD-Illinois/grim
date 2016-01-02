#include "reconstruction.hpp"

void reconstruction::reconstruct(const grid &prim,
                                 const int dir,
                                 grid &primLeft,
                                 grid &primRight
                                )
{
  switch (params::reconstruction)
  {
    case reconstructionOptions::MINMOD:

      reconstruction::reconstructMM(prim, dir, primLeft, primRight);

      break;

    case reconstructionOptions::WENO5:

      reconstruction::reconstructWENO5(prim, dir, primLeft, primRight);

      break;
  }

}
