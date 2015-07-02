#include "boundary.h"

void applyTileBoundaryConditions(const struct gridTile tile[ARRAY_ARGS 1],
                                 REAL primTile[ARRAY_ARGS TILE_SIZE(DOF)]
                                )
{
  /* 1) If PERIODIC, boundaries set by petsc during startFillingVecGhost
   * 2) if DIRICHLET, boundaries set by user in problem.c */
#if (   PHYSICAL_BOUNDARY_LEFT  != PERIODIC \
     && PHYSICAL_BOUNDARY_RIGHT != PERIODIC \
    )
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG,
                   -NG, TILE_SIZE_X2+NG,
                   -NG, TILE_SIZE_X1+NG
                  )
  {
    struct gridZone zone;
    setZone(iInTile, jInTile, kInTile, tile, &zone);

    /* LEFT EDGE */
    #if (PHYSICAL_BOUNDARY_LEFT != DIRICHLET)
    if (zone.iGlobal < 0)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_LEFT == OUTFLOW)
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(0, jInTile, kInTile, var)];
        #elif (PHYSICAL_BOUNDARY_LEFT == MIRROR)
          /* Mirror the left edge of the tile 
           * zone.iInTile goes from [-NG, 0) */
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(-iInTile-1, jInTile, kInTile, var)];
        #endif
      }
    }
    #endif /* Only compile when PHYSICAL_BOUNDARY_LEFT != DIRICHLET */
    /* END OF LEFT EDGE */

    /* RIGHT EDGE */
    #if (PHYSICAL_BOUNDARY_RIGHT != DIRICHLET)
    if (zone.iGlobal > N1-1)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_RIGHT == OUTFLOW)
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1, jInTile, kInTile, var)];
        #elif (PHYSICAL_BOUNDARY_RIGHT == MIRROR)
          /* zone.iInTile goes from [TILE_SIZE_X1, TILE_SIZE_X1+NG) 
           * iGhost goes from [0, NG) */
          int iGhost = iInTile - TILE_SIZE_X1;
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1-iGhost, jInTile, kInTile, var)];
        #endif
      }
    }
    #endif /* Only compile when PHYSICAL_BOUNDARY_RIGHT != DIRICHLET */
    /* END OF RIGHT EDGE */

  /* End of boundary conditions in the X1 direction */
  }
#endif /* Only compile when PHYSICAL_BOUNDARY != PERIODIC */

  /* It is necessary to make a second loop for the boundary conditions in X2 so
   * that the corner zones are also filled up */
#if (COMPUTE_DIM==2 || COMPUTE_DIM==3)
  #if (   PHYSICAL_BOUNDARY_TOP    != PERIODIC \
       && PHYSICAL_BOUNDARY_BOTTOM != PERIODIC \
      )
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG,
                   -NG, TILE_SIZE_X2+NG,
                   -NG, TILE_SIZE_X1+NG
                  )
  {
    struct gridZone zone;
    setZone(iInTile, jInTile, kInTile, tile, &zone);

    /* BOTTOM EDGE */
    #if (PHYSICAL_BOUNDARY_BOTTOM != DIRICHLET)
    if (zone.jGlobal < 0)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_BOTTOM == OUTFLOW)
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(iInTile, 0, kInTile, var)];
        #elif (PHYSICAL_BOUNDARY_BOTTOM == MIRROR)
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(iInTile, -jInTile-1, kInTile, var)];
        #endif
      }
    }
    #endif /* Only compile when PHYSICAL_BOUNDARY_BOTTOM != DIRICHLET */
    /* END OF BOTTOM EDGE */

    /* TOP EDGE */
    #if (PHYSICAL_BOUNDARY_TOP != DIRICHLET)
    if (zone.jGlobal > N2-1)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_TOP == OUTFLOW)
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(iInTile, TILE_SIZE_X2-1, kInTile, var)];
        #elif (PHYSICAL_BOUNDARY_TOP == MIRROR)
          int jGhost = jInTile - TILE_SIZE_X2;
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(iInTile, TILE_SIZE_X2-1-jGhost, kInTile, var)];
        #endif
      }
    }
    #endif /* Only compile when PHYSICAL_BOUNDARY_TOP != DIRICHLET */
    /* END OF TOP EDGE */

    /* Boundary conditions in the X2 direction */
  }
  #endif /* Only compile when PHYSICAL_BOUNDARY != PERIODIC */
#endif  /* Only compile when COMPUTE_DIM == 2 or 3 */

#if (COMPUTE_DIM==3)
  #if (   PHYSICAL_BOUNDARY_FRONT != PERIODIC \
       && PHYSICAL_BOUNDARY_BACK  != PERIODIC \
      )
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG,
                   -NG, TILE_SIZE_X2+NG,
                   -NG, TILE_SIZE_X1+NG
                  )
  {
    struct gridZone zone;
    setZone(iInTile, jInTile, kInTile, tile, &zone);

    /* FRONT */
    #if (PHYSICAL_BOUNDARY_FRONT != DIRICHLET)
    if (zone.kGlobal < 0)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_FRONT == OUTFLOW)
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(iInTile, jInTile, 0, var)];
        #elif (PHYSICAL_BOUNDARY_FRONT == MIRROR)
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(iInTile, jInTile, -kInTile-1, var)];
        #endif
      }
    }
    #endif /* Only compile when PHYSICAL_BOUNDARY_FRONT != DIRICHLET */
    /* END OF FRONT */

    /* BACK */
    #if (PHYSICAL_BOUNDARY_BACK != DIRICHLET)
    if (zone.kGlobal > N3-1)
    {
      for (int var=0; var<DOF; var++)
      {
        #if (PHYSICAL_BOUNDARY_BACK == OUTFLOW)
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(iInTile, jInTile, TILE_SIZE_X3-1, var)];
        #elif (PHYSICAL_BOUNDARY_BACK == MIRROR)
          int kGhost = kInTile - TILE_SIZE_X3;
          primTile[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE_MANUAL(iInTile, jInTile, TILE_SIZE_X3-1-kGhost, var)];
        #endif
      }
    }
    #endif /* Only compile when PHYSICAL_BOUNDARY_BACK != DIRICHLET */
    /* END OF BACK */

    /* Boundary conditions in the X3 direction */
  }
  #endif /* Only compile when PHYSICAL_BOUNDARY != PERIODIC */
#endif  /* Only compile when COMPUTE_DIM == 3 */

}
