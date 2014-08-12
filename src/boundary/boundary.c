#include "boundary.h"


void setZoneBoundaryFlags(struct gridZone zone[ARRAY_ARGS 1])
{
#if (COMPUTE_DIM==1 || COMPUTE_DIM==2)
  /* Zone is at the left end of a tile */
  if (zone->iTile == 0)
  {
    /* Is the zone at the left physical boundary? */
    if (zone->i == 0)
    {   
      /* If yes, then apply the physical boundary flag */
       zone->leftEdge = PHYSICAL_BOUNDARY_LEFT_EDGE; 
    }
    else
    {
      /* Zone is at the tile/MPI boundary but not at the physical boundary */
      zone->leftEdge = GET_DATA_FROM_LOCAL_ARRAY;
    }
  }

  /* Zone is at the right end of a tile */
  if (zone->iTile == TILE_SIZE_X1-1)
  {
    /* Is the zone at the right physical boundary? */
    if (zone->i == N1-1)
    {
        /* If yes, then apply the physical boundary condition flag */
        zone->rightEdge = PHYSICAL_BOUNDARY_RIGHT_EDGE;
    }
    else
    {   
        /* Zone is at the tile/MPI boundary but not at the physical boundary */
        zone->rightEdge = GET_DATA_FROM_LOCAL_ARRAY;
    }
  }
#endif /* Boundary flags in the X1 direction */

#if (COMPUTE_DIM==2)
    /* Zone is at the top end of a tile */
  if (zone->jTile == 0)
  {
    /* Is the zone at the top physical boundary? */
    if (zone->j == 0)
    {
      /* If yes, the apply the physical boundary flag */
      zone->topEdge = PHYSICAL_BOUNDARY_TOP_EDGE;
    }
    else
    {
      /* Zone is at the tile/MPI boundary but not at the physical boundary */
      zone->topEdge = GET_DATA_FROM_LOCAL_ARRAY;
    }
  }

  /* Zone is at the bottom end of a tile */
  if (zone->jTile == TILE_SIZE_X2-1)
  {
    /* Is the zone at the bottom physical boundary? */
    if (zone->j == N2-1)
    {
      /* If yes, then apply the physical boundary flag */
      zone->bottomEdge = PHYSICAL_BOUNDARY_BOTTOM_EDGE;
    }
    else
    {
      /* Zone is at the tile/MPI boundary but not at the physical boundary */
      zone->bottomEdge = GET_DATA_FROM_LOCAL_ARRAY;
    }
  }
#endif /* Boundary flags in the X2 direction */
}


void applyTileBoundaryConditions(const struct gridZone zone[ARRAY_ARGS 1],
                                 const REAL primLocal[ARRAY_ARGS TILE_SIZE],
                                 REAL tile[ARRAY_ARGS TILE_SIZE])
{
#if (COMPUTE_DIM==1 || COMPUTE_DIM==2)
    /* LEFT EDGE */
  if (zone->leftEdge == GET_DATA_FROM_LOCAL_ARRAY)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int iGhost=-NG; iGhost<0; iGhost++)
      {
        /* Get data from primLocal */
        tile[INDEX_TILE_MANUAL(iGhost, zone->jInTile, var)] =
        primLocal[INDEX_LOCAL_MANUAL(zone->i+iGhost, zone->j, zone, var)];
      }
    }
  }
  else if (zone->leftEdge == OUTFLOW)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int iGhost=-NG; iGhost<0; iGhost++)
      {
        tile[INDEX_TILE_MANUAL(iGhost, zone->jTile, var)] =
        tile[INDEX_TILE_MANUAL(0, zone->jTile, var)];
      }
    }
  }
  else if (zone->leftEdge == MIRROR)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int iGhost=-NG; iGhost<0; iGhost++)
      {
        /* Mirror the left edge of the tile */
        tile[INDEX_TILE_MANUAL(iGhost, zone->jTile, var)] = 
        tile[INDEX_TILE_MANUAL(-iGhost-1, zone->jTile, var)];
      }
    }
  }
  /* END OF LEFT EDGE */

    /* RIGHT EDGE */
  if (zone->rightEdge == GET_DATA_FROM_LOCAL_ARRAY)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int iGhost=0; iGhost<NG; iGhost++)
      {
        tile[INDEX_TILE_MANUAL(TILE_SIZE_X1+iGhost, zone->jInTile, var)] =
        primLocal[INDEX_LOCAL_MANUAL(zone->i+1+iGhost, zone->j, zone, var)];
      }
    }
  }
  else if (zone->rightEdge == OUTFLOW)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int iGhost=0; iGhost<NG; iGhost++)
      {
        tile[INDEX_TILE_MANUAL(TILE_SIZE_X1+iGhost, zone->jInTile, var)] =
        tile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1, zone->jInTile, var)];
      }
    }
  }
  else if (zone->rightEdge == MIRROR)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int iGhost=0; iGhost<NG; iGhost++)
      {
        tile[INDEX_TILE_MANUAL(TILE_SIZE_X1+iGhost, zone->jInTile, var)] =
        tile[INDEX_TILE_MANUAL(TILE_SIZE_X1-1-iGhost, zone->jInTile, var)];
      }
    }
  }
  /* END OF RIGHT EDGE */
#endif /* Boundary conditions in the X1 direction */

#if (COMPUTE_DIM==2)
    /* TOP EDGE */
  if (zone->topEdge == GET_DATA_FROM_LOCAL_ARRAY)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int jGhost=-NG; jGhost<0; jGhost++)
      {
        tile[INDEX_TILE_MANUAL(zone->iInTile, jGhost, var)] = 
        primLocal[INDEX_LOCAL_MANUAL(zone->i, zone->j+jGhost, zone, var)];
      }
    }
  }
  else if (zone->topEdge == OUTFLOW)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int jGhost=-NG; jGhost<0; jGhost++)
      {
        tile[INDEX_TILE_MANUAL(zone->iInTile, jGhost, var)] =
        tile[INDEX_TILE_MANUAL(zone->iInTile, 0, var)];
      }
    }
  }
  else if (zone->topEdge == MIRROR)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int jGhost=-NG; jGhost<0; jGhost++)
      {
        tile[INDEX_TILE_MANUAL(zone->iInTile, jGhost, var)] = 
        tile[INDEX_TILE_MANUAL(zone->iInTile, -jGhost-1, var)];
      }
    } 
  }
  /* END OF TOP EDGE */

  /* BOTTOM EDGE */
  if (zone->bottomEdge == GET_DATA_FROM_LOCAL_ARRAY)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int jGhost=0; jGhost<NG; jGhost++)
      {
        tile[INDEX_TILE_MANUAL(zone->iInTile, TILE_SIZE_X2+jGhost, var)] =
        primLocal[INDEX_LOCAL_MANUAL(zone->i, zone->j+1+jGhost, zone, var)];
      }
    }
  }
  else if (zone->bottomEdge == OUTFLOW)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int jGhost=0; jGhost<NG; jGhost++)
      {
        tile[INDEX_TILE_MANUAL(zone->iInTile, TILE_SIZE_X2+jGhost, var)] = 
        tile[INDEX_TILE_MANUAL(zone->iInTile, TILE_SIZE_X2-1, var)];
      }
    }
  }
  else if (zone->bottomEdge == MIRROR)
  {
    for (int var=0; var<DOF; var++)
    {
      for (int jGhost=0; jGhost<NG; jGhost++)
      {
        tile[INDEX_TILE_MANUAL(zone->iInTile, TILE_SIZE_X2+jGhost, var)] = 
        tile[INDEX_TILE_MANUAL(zone->iInTile, TILE_SIZE_X2-1-jGhost, var)];
      }
    }
  }
  /* END OF BOTTOM EDGE */
#endif /* Boundary conditions in the X2 direction */
}
