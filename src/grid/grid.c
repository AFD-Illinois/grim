#include "grid.h"

/* Get the coordinates corresponding to a location inside a grid zone.
 *         
 * 2D:
 *           +------------+
 *           |            |
 *           |            |
 * (i,j+0.5) F1    C      |
 *           |            |
 *           |            |
 *           +-----F2-----+
 *         (i,j)  (i+0.5,j)
 *
 *         X2
 *         ^
 *         |
 *         |
 *         +---> X1
 *      C   --  Center (CENTER)
 *      F1  --  Face in the X1 direction (FACE_X1)
 *      F2  --  Face in the X2 direction (FACE_X2)
 *
 * 3D:
 *
 *
 *           +------+  
 *          /|     /|  
 *         +-+----+ |  
 *         | |    | |  
 *         | +----+-+  
 *         |/     |/   
 *         +------+    
 *      (i,j,k)   
 *
 *      X2
 *      ^
 *      |
 *      |
 *      +---> X1
 *     /
 *    X3
 *      Center (CENTER)                    (i+0.5, j+0.5, k+0.5)
 *      Face in the X1 direction (FACE_X1) (i, j+0.5, k+0.5)
 *      Face in the X2 direction (FACE_X2) (i+0.5, j, k+0.5)
 *      Face in the X3 direction (FACE_X3) (i+0.5, j+0.5, k) 
 *      Corner (CORNER)                    (i,j,k)
 *
 * @param input: zone, Zone in which the coordinates are needed
 * @param input: location, Location inside the zone where the coordinates are needed
 * @param output: X, The coordinates X^mu, in which the computation is
 *                performed.
*/
void getXCoords(const struct gridZone zone[ARRAY_ARGS 1],
                const int location,
                REAL X[ARRAY_ARGS NDIM])
{
  X[0] = 0.; /* Don't care about the time coordinate. Just set to 0.*/

  switch (location)
  {
    case CORNER:            /* (i,j,k) */
      
      X[1] = X1_A + zone->iGlobal*zone->dX[1];
      X[2] = X2_A + zone->jGlobal*zone->dX[2];
      X[3] = X3_A + zone->kGlobal*zone->dX[3];

      break;

    case CENTER:            /* (i+0.5,j+0.5, k+0.5) */

      X[1] = X1_A + (zone->iGlobal + 0.5)*zone->dX[1];
      X[2] = X2_A + (zone->jGlobal + 0.5)*zone->dX[2];
      X[3] = X3_A + (zone->kGlobal + 0.5)*zone->dX[3];

      break;

    case FACE_X1:           /* (i,j+0.5,k+0.5) */

      X[1] = X1_A + (zone->iGlobal      )*zone->dX[1];
      X[2] = X2_A + (zone->jGlobal + 0.5)*zone->dX[2];
      X[3] = X3_A + (zone->kGlobal + 0.5)*zone->dX[3];

      break;

    case FACE_X1_PLUS_ONE:  /* (i+1,j+0.5,k+0.5) */

      X[1] = X1_A + (zone->iGlobal + 1  )*zone->dX[1];
      X[2] = X2_A + (zone->jGlobal + 0.5)*zone->dX[2];
      X[3] = X3_A + (zone->kGlobal + 0.5)*zone->dX[3];

      break;

    case FACE_X2:           /* (i+0.5,j,k+0.5) */

      X[1] = X1_A + (zone->iGlobal + 0.5)*zone->dX[1];
      X[2] = X2_A + (zone->jGlobal      )*zone->dX[2];
      X[3] = X3_A + (zone->kGlobal + 0.5)*zone->dX[3];

      break;

    case FACE_X2_PLUS_ONE:  /* (i+0.5,j+1,k+0.5) */

      X[1] = X1_A + (zone->iGlobal + 0.5)*zone->dX[1];
      X[2] = X2_A + (zone->jGlobal + 1  )*zone->dX[2];
      X[3] = X3_A + (zone->kGlobal + 0.5)*zone->dX[3];

      break;

    case FACE_X3:           /* (i+0.5,j+0.5,k) */

      X[1] = X1_A + (zone->iGlobal + 0.5)*zone->dX[1];
      X[2] = X2_A + (zone->jGlobal + 0.5)*zone->dX[2];
      X[3] = X3_A + (zone->kGlobal      )*zone->dX[3];

      break;

    case FACE_X3_PLUS_ONE:  /* (i+0.5,j+0.5,k+1) */

      X[1] = X1_A + (zone->iGlobal + 0.5)*zone->dX[1];
      X[2] = X2_A + (zone->jGlobal + 0.5)*zone->dX[2];
      X[3] = X3_A + (zone->kGlobal + 1  )*zone->dX[3];

      break;

    default:
      
      #if (DEBUG)

      #endif
      break;
  }
}

void initGridData(const int numVar, const int numGhost,
                  struct gridData grid[ARRAY_ARGS 1]
                 )
{
  /* Periodic boundary conditions handled by Petsc since it is a global boundary
   * condition. Here we check for the boundary at the left edge. Obviously the
   * boundary at the right edge also must be PERIODIC if left edge is PERIODIC */
  DMBoundaryType boundaryX1 = DM_BOUNDARY_GHOSTED;
  DMBoundaryType boundaryX2 = DM_BOUNDARY_GHOSTED;
  DMBoundaryType boundaryX3 = DM_BOUNDARY_GHOSTED;

  #if (PHYSICAL_BOUNDARY_LEFT == PERIODIC)
    boundaryX1 = DM_BOUNDARY_PERIODIC;
  #endif

  #if (PHYSICAL_BOUNDARY_TOP == PERIODIC)
    boundaryX2 = DM_BOUNDARY_PERIODIC;
  #endif

  #if (PHYSICAL_BOUNDARY_FRONT == PERIODIC)
    boundaryX3 = DM_BOUNDARY_PERIODIC;
  #endif

  #if (COMPUTE_DIM==1)

    DMDACreate1d(PETSC_COMM_WORLD, boundaryX1, N1,
                 numVar, 0, NULL,
                 &grid->dm
                );

    DMDACreate1d(PETSC_COMM_WORLD, boundaryX1, N1,
                 numVar, numGhost, NULL,
                 &grid->dmGhost
                );
  
  #elif (COMPUTE_DIM==2)

    DMDACreate2d(PETSC_COMM_WORLD, 
                 boundaryX1, boundaryX2,
                 DMDA_STENCIL_BOX,
                 N1, N2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 numVar, 0,
                 PETSC_NULL, PETSC_NULL,
                 &grid->dm
                );

    DMDACreate2d(PETSC_COMM_WORLD, 
                 boundaryX1, boundaryX2,
                 DMDA_STENCIL_BOX,
                 N1, N2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 numVar, numGhost,
                 PETSC_NULL, PETSC_NULL,
                 &grid->dmGhost
                );

  #elif (COMPUTE_DIM==3)

    DMDACreate3d(PETSC_COMM_WORLD, 
                 boundaryX1, boundaryX2, boundaryX3,
                 DMDA_STENCIL_BOX,
                 N1, N2, N3,
                 PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE
                 numVar, 0, 
                 PETSC_NULL, PETSC_NULL, PETSC_NULL,
                 &grid->dm
                );

    DMDACreate3d(PETSC_COMM_WORLD, 
                 boundaryX1, boundaryX2, boundaryX3,
                 DMDA_STENCIL_BOX,
                 N1, N2, N3,
                 PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE
                 numVar, numGhost, 
                 PETSC_NULL, PETSC_NULL, PETSC_NULL,
                 &grid->dmGhost
                );

  #endif /* Choose dim and create dm */

  /* Get local indices for looping over the array owned by the current proc */
  DMDAGetCorners(grid->dm,
                 &grid->iLocalStart, &grid->jLocalStart, &grid->kLocalStart,
                 &grid->iLocalSize,  &grid->jLocalSize,  &grid->kLocalSize
                );

  grid->numTilesX1 = grid->iLocalSize/TILE_SIZE_X1;
  grid->numTilesX2 = grid->jLocalSize/TILE_SIZE_X2;

  DMCreateGlobalVector(grid->dm, &grid->vec);
  DMCreateLocalVector(grid->dmGhost, &grid->vecGhost);

  VecSet(grid->vec, 0.);
  VecSet(grid->vecGhost, 0.);
}

void destroyGridData(struct gridData grid[ARRAY_ARGS 1])
{
  VecDestroy(&grid->vec);
  VecDestroy(&grid->vecGhost);

  DMDestroy(&grid->dm);
  DMDestroy(&grid->dmGhost);
}

void setPointerToVec(struct gridData grid[ARRAY_ARGS 1])
{
  DMDAVecGetArrayDOF(grid->dm, grid->vec, &grid->ptr);
}

void restorePointerToVec(struct gridData grid[ARRAY_ARGS 1])
{
  DMDAVecRestoreArrayDOF(grid->dm, grid->vec, &grid->ptr);
}

void setPointerToVecGhost(struct gridData grid[ARRAY_ARGS 1])
{
  DMDAVecGetArrayDOF(grid->dm, grid->vec, &grid->ptrGhost);
}

void restorePointerToVecGhost(struct gridData grid[ARRAY_ARGS 1])
{
  DMDAVecRestoreArrayDOF(grid->dm, grid->vec, &grid->ptrGhost);
}

void setPointerToExternalVec(Vec externalVec,
                             struct gridData grid[ARRAY_ARGS 1]
                            )
{
  DMDAVecGetArrayDOF(grid->dm, externalVec, &grid->ptr);
}

void restorePointerToExternalVec(Vec externalVec,
                                 struct gridData grid[ARRAY_ARGS 1]
                                )
{
  DMDAVecRestoreArrayDOF(grid->dm, externalVec, &grid->ptr);
}

void startFillingVecGhost(struct gridData grid[ARRAY_ARGS 1])
{
  DMGlobalToLocalBegin(grid->dm, grid->vec, INSERT_VALUES, grid->vecGhost); 
}

void finishFillingVecGhost(struct gridData grid[ARRAY_ARGS 1])
{
  DMGlobalToLocalEnd(grid->dm, grid->vec, INSERT_VALUES, grid->vecGhost); 
}

void startFillingVecGhostWithExternalVec(Vec externalVec,
                                         struct gridData grid[ARRAY_ARGS 1]
                                        )
{
  DMGlobalToLocalBegin(grid->dm, externalVec, INSERT_VALUES, grid->vecGhost); 
}

void finishFillingVecGhostWithExternalVec(Vec externalVec,
                                          struct gridData grid[ARRAY_ARGS 1]
                                         )
{
  DMGlobalToLocalEnd(grid->dm, externalVec, INSERT_VALUES, grid->vecGhost); 
}

/* Function to set a {\tt gridZone}.
 *
 * Note that the macro LOOP_OVER_TILES defines the variables iTile, jTile and 
 * kGlobal.
 *
 * @param: Input: iInTile, X1 index inside a tile.
 * @param: Input: jInTile, X2 index inside a tile.
 * @param: Input: tile, a {\tt gridTile} struct.
 * @param: Output: zone, the {\tt gridZone} with its parameters set.
 */
void setZone(const int iInTile, const int jInTile, const int kInTile,
             const struct gridTile tile[ARRAY_ARGS 1],
             struct gridZone zone[ARRAY_ARGS 1]
            )
{
  zone->iInTile = iInTile;
  zone->jInTile = jInTile;
  zone->kInTile = kInTile;

  zone->iGlobal = tile->iLocalStart + iInTile + tile->iTile*TILE_SIZE_X1;
  zone->jGlobal = tile->jLocalStart + jInTile + tile->jTile*TILE_SIZE_X2;
  zone->kGlobal = kInTile + tile->kGlobal; /* kInTile = 0 == kGlobal */
  
  zone->dX[0] = 0.;
  zone->dX[1] = (X1_B - X1_A)/((REAL)N1);
  #if (COMPUTE_DIM==1)
    zone->dX[2] = 0.;
    zone->dX[3] = 0.;
  #elif (COMPUTE_DIM==2)
    zone->dX[2] = (X2_B - X2_A)/((REAL)N2);
    zone->dX[3] = 0.;
  #elif (COMPUTE_DIM==3)
    zone->dX[2] = (X2_B - X2_A)/((REAL)N2);
    zone->dX[3] = (X3_B - X3_A)/((REAL)N3);
  #endif
}

void setTile(const int iTile, const int jTile, const int kGlobal,
             const struct gridData grid[ARRAY_ARGS 1],
             struct gridTile tile[ARRAY_ARGS 1]
            )
{
  tile->iTile = iTile;
  tile->jTile = jTile;
  tile->kGlobal = kGlobal;

  tile->iLocalStart = grid->iLocalStart;
  tile->jLocalStart = grid->jLocalStart;
  tile->kLocalStart = grid->kLocalStart;
}

void loadPrimTile(const struct gridData grid[ARRAY_ARGS 1],
                  const struct gridTile tile[ARRAY_ARGS 1],
                  REAL primTile[ARRAY_ARGS TILE_SIZE(DOF)]
                 )
{
  if (tile->kGlobal == tile->kLocalStart)
  {
    LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG,
                     -NG, TILE_SIZE_X2+NG,
                     -NG, TILE_SIZE_X3+NG
                    )
    {
      struct gridZone zone;
      setZone(iInTile, jInTile, kInTile, tile, &zone);

      REAL primVars[DOF];
      memcpy(primVars, &INDEX_GRID_GHOST(grid, &zone, 0), 
             sizeof(REAL[DOF])
            );

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] = primVars[var];
      }
    }
  }
  else
  {
    struct gridZone zoneAtTileEdge;
    setZone(0, 0, 0, tile, &zoneAtTileEdge);

    memmove(&primTile[INDEX_TILE_OFFSET(0, 0, NG,   &zoneAtTileEdge, 0)],
            &primTile[INDEX_TILE_OFFSET(0, 0, NG-1, &zoneAtTileEdge, 0)],
            sizeof(REAL[(TILE_SIZE_X1+2*NG)*(TILE_SIZE_X2+2*NG)*(2*NG)*DOF])
           );

    LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG,
                     -NG, TILE_SIZE_X2+NG,
                      NG, TILE_SIZE_X3+NG
                    )
    {
      struct gridZone zone;
      setZone(iInTile, jInTile, kInTile, tile, &zone);

      REAL primVars[DOF];
      memcpy(primVars, &INDEX_GRID_GHOST(grid, &zone, 0), 
             sizeof(REAL[DOF])
            );

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] = primVars[var];
      }
    }
  }
}
