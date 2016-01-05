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

array reconstruction::slope(const int dir,const double dX,
			    const array& in)
{
  switch (params::reconstruction)
    {
    case reconstructionOptions::MINMOD:
      return reconstruction::slopeMM(dir,dX,in);
    case reconstructionOptions::WENO5:
      return reconstruction::slopeWENO5(dir,dX,in);
   }
}
