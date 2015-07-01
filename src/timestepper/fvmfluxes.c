#include "timestepper.h"

void computeFluxesOverTile
  (const REAL primTile[ARRAY_ARGS TILE_SIZE(DOF)],
   const struct gridTile tile[ARRAY_ARGS 1],
   REAL fluxesTile[COMPUTE_DIM][ARRAY_ARGS TILE_SIZE(DOF)],
   struct gridData dtGrid[ARRAY_ARGS 1]
  )
{
  REAL primTileLeft[TILE_SIZE(DOF)], primTileRight[TILE_SIZE(DOF)];
  
  REAL dt[COMPUTE_DIM];
  dt[0]=1e10, dt[1]=1e10, dt[2]=1e10;

  int indexFace[4] = {0, FACE_X1, FACE_X2, FACE_X3};
  int iOffset[4]   = {0, -1, 0, 0};     
  int jOffset[4]   = {0, 0, -1, 0};     
  int kOffset[4]   = {0, 0, 0, -1};     

  for (int dir=1; dir <= COMPUTE_DIM; dir++)
  {
    reconstruct(primTile, tile, dir,
                primVarsLeft, primVarsRight
               );

    LOOP_INSIDE_TILE(-1, TILE_SIZE_X1+1,
                     -1, TILE_SIZE_X2+1,
                     -1, TILE_SIZE_X3+1
                    )
    {
      struct gridZone zone;
      setGridZone(iInTile, jInTile, kInTile, tile, &zone);

      REAL primVarsLeft[DOF]      , primVarsRight[DOF];
      REAL fluxesLeft[DOF]        , fluxesRight[DOF];
      REAL conservedVarsLeft[DOF] , conservedVarsRight[DOF];

      memcpy
        (primVarsLeft, 
         primTileRight[INDEX_TILE_OFFSET(iOffset[dir], 
                                         jOffset[dir],
                                         kOffset[dir],
                                         &zone, 0
                                        )
                      ],
         sizeof(REAL[DOF])
        );

      REAL XCoords[NDIM];
      getXCoords(&zone, indexFace[dir], XCoords);
      struct geometry geom; setGeometry(XCoords, &geom);
      struct fluidElement elem;
      setFluidElement(primVarsleft, &geom, &elem);

      computeFluxes(&elem, &geom, dir, fluxesLeft);
      computeFluxes(&elem, &geom, 0  , conservedVarsLeft);

      memcpy
        (primVarsRight, 
         primTileLeft[INDEX_TILE(&zone, 0)],
         sizeof(REAL[DOF])
        );

      setFluidElement(primVarsRight, &geom, &elem);

      computeFluxes(&elem, &geom, dir, fluxesRight);
      computeFluxes(&elem, &geom, 0  , conservedVarsRight);

      REAL waveSpeed =
        riemannSolver(fluxesLeft        , fluxesRight,
                      conservedVarsLeft , conservedVarsRight,
                      primVarsLeft      , primVarsRight,
                      &geom, dir, &fluxesTile[dir-1][INDEX_TILE(&zone, 0)]
                     ); 

      /* Compute the dt in each zone */
      zone.dX[0] = COURANT * zone.dX[dir] / waveSpeed;

      if (zone.dX[0] < dt[dir-1])
      {
        dt[dir-1] = zone.dX[0];
      }
    }

    struct gridZone zoneAtTileEdge;
    setZone(0, 0, 0, tile, &zoneAtTileEdge);
    memset(INDEX_GRID(dtGrid, &zoneAtTileEdge, dir-1),
           dt[dir-1], sizeof(REAL[TILE_SIZE(1)])
          );
  }
        
  /* Flux CT */
  fluxCT(fluxesTile);

}
