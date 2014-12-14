#include "boundary.h"

void applyTileBoundaryConditions(const int iTile, const int jTile,
                                 const int X1Start, const int X2Start,
                                 const int X1Size, const int X2Size,
                                 REAL tile[ARRAY_ARGS TILE_SIZE])
{
  /* First loop over all ghost zones, including corners.  This will set the
   * values of all ghost zones except the corners of the global domain. */
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

    /* LEFT EDGE */
    if (zone.i < 0)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_LEFT_EDGE == OUTFLOW)
          tile[INDEX_TILE(&zone, var)] =
          tile[INDEX_TILE_MANUAL(0, zone.jInTile, var)];
        #elif (PHYSICAL_BOUNDARY_LEFT_EDGE == MIRROR)
          /* Mirror the left edge of the tile 
           * zone.iInTile goes from [-NG, 0) */
          tile[INDEX_TILE(&zone, var)] =
          tile[INDEX_TILE_MANUAL(-zone.iInTile-1, zone.jInTile, var)];
        #endif
      }
    }
    /* END OF LEFT EDGE */


    /* RIGHT EDGE */
    if (zone.i > N1-1)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_RIGHT_EDGE == OUTFLOW)
          tile[INDEX_TILE(&zone, var)] =
          tile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1, zone.jInTile, var)];
        #elif (PHYSICAL_BOUNDARY_RIGHT_EDGE == MIRROR)
          /* zone.iInTile goes from [TILE_SIZE_X1, TILE_SIZE_X1+NG) 
           * iGhost goes from [0, NG) */
          int iGhost = zone.iInTile - TILE_SIZE_X1;
          tile[INDEX_TILE(&zone, var)] =
          tile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1-iGhost, zone.jInTile, var)];
        #endif
      }
    }
    /* END OF RIGHT EDGE */

  /* End of boundary conditions in the X1 direction */

  }

#if (COMPUTE_DIM==2)
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

    /* BOTTOM EDGE */
    if (zone.j < 0)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_BOTTOM_EDGE == OUTFLOW)
          tile[INDEX_TILE(&zone, var)] =
          tile[INDEX_TILE_MANUAL(zone.iInTile, 0, var)];
        #elif (PHYSICAL_BOUNDARY_BOTTOM_EDGE == MIRROR)
          tile[INDEX_TILE(&zone, var)] =
          tile[INDEX_TILE_MANUAL(zone.iInTile, -zone.jInTile-1, var)];
        #endif
      }
    }
    /* END OF BOTTOM EDGE */

    /* TOP EDGE */
    if (zone.j > N2-1)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_TOP_EDGE == OUTFLOW)
          tile[INDEX_TILE(&zone, var)] =
          tile[INDEX_TILE_MANUAL(zone.iInTile, TILE_SIZE_X2-1, var)];
        #elif (PHYSICAL_BOUNDARY_TOP_EDGE == MIRROR)
          int jGhost = zone.jInTile - TILE_SIZE_X2;
          tile[INDEX_TILE(&zone, var)] =
          tile[INDEX_TILE_MANUAL(zone.iInTile, TILE_SIZE_X2-1-jGhost, var)];
        #endif
      }
    }
    /* END OF TOP EDGE */

    /* Boundary conditions in the X2 direction */
  }
#endif

//#if (COMPUTE_DIM==2)
//  /* Now loop over the corner zones of the global domain */
//  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
//  {
//    struct gridZone zone;
//    setGridZone(iTile, jTile,
//                iInTile, jInTile,
//                X1Start, X2Start, 
//                X1Size, X2Size, 
//                &zone);
//
//    /* LEFT EDGE */
//    if (zone.i < 0 && (zone.j < 0 || zone.j > N2-1) )
//    {
//      for (int var=0; var<DOF; var++)
//      {
//        #if (PHYSICAL_BOUNDARY_LEFT_EDGE == OUTFLOW)
//          tile[INDEX_TILE(&zone, var)] =
//          tile[INDEX_TILE_MANUAL(0, zone.jInTile, var)];
//        #elif (PHYSICAL_BOUNDARY_LEFT_EDGE == MIRROR)
//          /* Mirror the left edge of the tile 
//           * zone.iInTile goes from [-NG, 0) */
//          tile[INDEX_TILE(&zone, var)] =
//          tile[INDEX_TILE_MANUAL(-zone.iInTile-1, zone.jInTile, var)];
//        #endif
//      }
//    }
//    /* END OF LEFT EDGE */
//
//
//    /* RIGHT EDGE */
//    if (zone.i > N1-1 && (zone.j < 0 || zone.j > N2-1) )
//    {
//      for (int var=0; var<DOF; var++)
//      {
//        #if (PHYSICAL_BOUNDARY_RIGHT_EDGE == OUTFLOW)
//          tile[INDEX_TILE(&zone, var)] =
//          tile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1, zone.jInTile, var)];
//        #elif (PHYSICAL_BOUNDARY_RIGHT_EDGE == MIRROR)
//          /* zone.iInTile goes from [TILE_SIZE_X1, TILE_SIZE_X1+NG) 
//           * iGhost goes from [0, NG) */
//          int iGhost = zone.iInTile - TILE_SIZE_X1;
//          tile[INDEX_TILE(&zone, var)] =
//          tile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1-iGhost, zone.jInTile, var)];
//        #endif
//      }
//    }
//    /* END OF RIGHT EDGE */
//
//    /* End of boundary conditions of the corner zones */
//  }
//#endif
}
