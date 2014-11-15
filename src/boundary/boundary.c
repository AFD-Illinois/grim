#include "boundary.h"

void applyTileBoundaryConditions(const int iTile, const int jTile,
                                 const int X1Start, const int X2Start,
                                 const int X1Size, const int X2Size,
                                 ARRAY(primLocal),
                                 REAL tile[ARRAY_ARGS TILE_SIZE])
{
  LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

    /* LEFT EDGE */
    if (zone.i == 0)
    {
      for (int var=0; var<DOF; var++)
      {
        for (int iGhost=-NG; iGhost<0; iGhost++)
        {
          #if (PHYSICAL_BOUNDARY_LEFT_EDGE == OUTFLOW)
            tile[INDEX_TILE_MANUAL(iGhost, zone.jTile, var)] =
            tile[INDEX_TILE_MANUAL(0, zone.jTile, var)];
          #elif (PHYSICAL_BOUNDARY_LEFT_EDGE == MIRROR)
            /* Mirror the left edge of the tile */
            tile[INDEX_TILE_MANUAL(iGhost, zone.jTile, var)] = 
            tile[INDEX_TILE_MANUAL(-iGhost-1, zone.jTile, var)];
          #endif
        }
      }
    }
    /* END OF LEFT EDGE */


    /* RIGHT EDGE */
    if (zone.i == N1-1)
    {
      for (int var=0; var<DOF; var++)
      {
        for (int iGhost=0; iGhost<NG; iGhost++)
        {
          #if (PHYSICAL_BOUNDARY_RIGHT_EDGE == OUTFLOW)
            tile[INDEX_TILE_MANUAL(TILE_SIZE_X1+iGhost, zone.jInTile, var)] =
            tile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1, zone.jInTile, var)];
          #elif (PHYSICAL_BOUNDARY_RIGHT_EDGE == MIRROR)
            tile[INDEX_TILE_MANUAL(TILE_SIZE_X1+iGhost, zone.jInTile, var)] =
            tile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1-iGhost, zone.jInTile, var)];
          #endif
        }
      }
    }
    /* END OF RIGHT EDGE */

  /* End of boundary conditions in the X1 direction */

#if (COMPUTE_DIM==2)
    /* BOTTOM EDGE */
    if (zone.j == 0)
    {
      for (int var=0; var<DOF; var++)
      {
        for (int jGhost=-NG; jGhost<0; jGhost++)
        {
          #if (PHYSICAL_BOUNDARY_BOTTOM_EDGE == OUTFLOW)
            tile[INDEX_TILE_MANUAL(zone.iInTile, jGhost, var)] =
            tile[INDEX_TILE_MANUAL(zone.iInTile, 0, var)];
          #elif (PHYSICAL_BOUNDARY_BOTTOM_EDGE == MIRROR)
            tile[INDEX_TILE_MANUAL(zone.iInTile, jGhost, var)] = 
            tile[INDEX_TILE_MANUAL(zone.iInTile, -jGhost-1, var)];
          #endif
        }
      }
    }
    /* END OF BOTTOM EDGE */

    /* TOP EDGE */
    if (zone.j == N2-1)
    {
      for (int var=0; var<DOF; var++)
      {
        for (int jGhost=0; jGhost<NG; jGhost++)
        {
          #if (PHYSICAL_BOUNDARY_TOP_EDGE == OUTFLOW)
            tile[INDEX_TILE_MANUAL(zone.iInTile, TILE_SIZE_X2+jGhost, var)] = 
            tile[INDEX_TILE_MANUAL(zone.iInTile, TILE_SIZE_X2-1, var)];
          #elif (PHYSICAL_BOUNDARY_TOP_EDGE == MIRROR)
            tile[INDEX_TILE_MANUAL(zone.iInTile, TILE_SIZE_X2+jGhost, var)] = 
            tile[INDEX_TILE_MANUAL(zone.iInTile, TILE_SIZE_X2-1-jGhost, var)];
          #endif
        }
      }
    }
    /* END OF TOP EDGE */

#endif /* Boundary conditions in the X2 direction */
  }

}
