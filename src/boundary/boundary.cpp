#include "boundary.hpp"

namespace boundaries
{
  void applyBoundaryConditions(const int boundaryLeft, const int boundaryRight,
                               const int boundaryTop, const int boundaryBottom,
                               const int boundaryFront, const int boundaryBack,
                               grid &prim
                              )
  {
    int numGhost = prim.numGhost;
    
    /* Domain of X1 boundaries */
    af::seq domainX1LeftBoundary(0, numGhost-1);
    af::seq domainX1RightBoundary(prim.N1Local+numGhost, 
                                  prim.N1Local+2*numGhost-1
                                 );
    
    /* Domain of X2 boundaries */
    af::seq domainX2BottomBoundary(0, numGhost-1);
    af::seq domainX2TopBoundary(prim.N2Local+numGhost,
                                prim.N2Local+2*numGhost-1
                               );

    /* Domain of X3 boundaries */
    af::seq domainX3BackBoundary(0, numGhost-1);
    af::seq domainX3FrontBoundary(prim.N3Local+numGhost,
                                  prim.N3Local+2*numGhost-1
                                 );


    /* X1 boundary conditions */
    if (prim.iLocalStart == 0)
    {
      switch (boundaryLeft)
      {
        case (boundaries::MIRROR):
        {
          af::seq domainX1EdgeReversed(2*numGhost-1, numGhost, -1);

          for (int var=0; var<prim.numVars; var++)
          {
            prim.vars[var](domainX1LeftBoundary, span, span) =
              prim.vars[var](domainX1EdgeReversed, span, span);
          }

          break;
        }

        case (boundaries::OUTFLOW):
        {
          for (int var=0; var<prim.numVars; var++)
          {
            gfor (af::seq i, 0, numGhost-1)
            {
              prim.vars[var](i, span, span) =
                prim.vars[var](numGhost, span, span);
            }
          }

          break;
        }
      }
    }

    if (prim.iLocalEnd == prim.N1)
    {
      switch (boundaryRight)
      {
        case (boundaries::MIRROR):
        {
          af::seq domainX1EdgeReversed(prim.N1Local+numGhost-1, prim.N1Local, -1);

          for (int var=0; var<prim.numVars; var++)
          {
            prim.vars[var](domainX1RightBoundary, span, span) =
              prim.vars[var](domainX1EdgeReversed, span, span);
          }

          break;
        }

        case (boundaries::OUTFLOW):
        {
          for (int var=0; var<prim.numVars; var++)
          {
            gfor (af::seq i, prim.N1Local+numGhost, 
                             prim.N1Local+2*numGhost-1
                 )
            {
              prim.vars[var](i, span, span) =
                prim.vars[var](prim.N1Local+numGhost-1, span, span);
            }
          }

          break;
        }
      }
    }

    /* X2 boundary conditions */
    if (prim.jLocalStart == 0 && prim.dim >= 2)
    {
      switch (boundaryBottom)
      {
        case (boundaries::MIRROR):
        {
          af::seq domainX2EdgeReversed(2*numGhost-1, numGhost, -1);

          for (int var=0; var<prim.numVars; var++)
          {
            prim.vars[var](span, domainX2BottomBoundary, span)
              = prim.vars[var](span, domainX2EdgeReversed, span);
          }

          break;
        }

        case (boundaries::OUTFLOW):
        {
          for (int var=0; var<prim.numVars; var++)
          {
            gfor (af::seq j, 0, numGhost-1)
            {
              prim.vars[var](span, j, span)
                = prim.vars[var](span, numGhost, span);
            }
          }

          break;
        }
      }
    }

    if (prim.jLocalEnd == prim.N2 && prim.dim >= 2)
    {
      switch (boundaryTop)
      {
        case (boundaries::MIRROR):
        {
          af::seq domainX2EdgeReversed(prim.N2Local+numGhost-1, prim.N2Local, -1);

          for (int var=0; var<prim.numVars; var++)
          {
            prim.vars[var](span, domainX2TopBoundary, span) =
              prim.vars[var](span, domainX2EdgeReversed, span);
          
          }

          break;
        }

        case (boundaries::OUTFLOW):
        {
          for (int var=0; var<prim.numVars; var++)
          {
            gfor (af::seq j, prim.N2Local+numGhost, 
                             prim.N2Local+2*numGhost-1
                 )
            {
              prim.vars[var](span, j, span) =
                prim.vars[var](span, prim.N2Local+numGhost-1, span);
            }
          }

          break;
        }
      }
    }

    /* X3 boundary conditions */
    if (prim.kLocalStart == 0 && prim.dim >= 3)
    {
      switch (boundaryBack)
      {
        case (boundaries::MIRROR):
        {
          af::seq domainX3EdgeReversed(2*numGhost-1, numGhost, -1);

          for (int var=0; var<prim.numVars; var++)
          {
            prim.vars[var](span, span, domainX3BackBoundary)
              = prim.vars[var](span, span, domainX3EdgeReversed);
          }

          break;
        }

        case (boundaries::OUTFLOW):
        {
          for (int var=0; var<prim.numVars; var++)
          {
            gfor (af::seq k, 0, numGhost-1)
            {
              prim.vars[var](span, span, k)
                = prim.vars[var](span, span, numGhost);
            }
          }

          break;
        }
      }
    }

    if (prim.kLocalEnd == prim.N3 && prim.dim >= 3)
    {
      switch (boundaryFront)
      {
        case (boundaries::MIRROR):
        {
          af::seq domainX3EdgeReversed(prim.N3Local+numGhost-1, prim.N3Local, -1);

          for (int var=0; var<prim.numVars; var++)
          {
            prim.vars[var](span, span, domainX3FrontBoundary) =
              prim.vars[var](span, span, domainX3EdgeReversed);
          
          }

          break;
        }

        case (boundaries::OUTFLOW):
        {
          for (int var=0; var<prim.numVars; var++)
          {
            gfor (af::seq k, prim.N3Local+numGhost, 
                             prim.N3Local+2*numGhost-1
                 )
            {
              prim.vars[var](span, span, k) =
                prim.vars[var](span, span, prim.N3Local+numGhost-1);
            }
          }

          break;
        }
      }
    }

  } /* End of applyBoundaryConditions() */

} /* namespace boundaries */
