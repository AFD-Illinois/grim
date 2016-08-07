#include "reconstruction.hpp"

void reconstruction::reconstruct(const grid &prim,
                                 const int dir,
                                 grid &primLeft,
                                 grid &primRight,
                                 int &numReads,
                                 int &numWrites
                                )
{
  switch (params::reconstruction)
  {
    case reconstructionOptions::MINMOD:

      reconstruction::reconstructMM(prim, dir, primLeft, primRight,
                                    numReads, numWrites
                                   );

      break;

    case reconstructionOptions::WENO5:

      reconstruction::reconstructWENO5(prim, dir, primLeft, primRight,
                                       numReads, numWrites
                                      );

      break;

  case reconstructionOptions::PPM:

      reconstruction::reconstructPPM(prim, dir, primLeft, primRight,
                                       numReads, numWrites
                                      );

      break;
  }

}

array reconstruction::slope(const int dir, const double dX,
			                      const array& in,
                            int &numReads,
                            int &numWrites
                           )
{
  switch (params::reconstruction)
  {
    case reconstructionOptions::MINMOD:

      return reconstruction::slopeMM(dir,dX,in, numReads, numWrites);

    case reconstructionOptions::WENO5:

      return reconstruction::slopeWENO5(dir,dX,in, numReads, numWrites);

    case reconstructionOptions::PPM:

      return reconstruction::slopePPM(dir,dX,in, numReads, numWrites);
   }
}
