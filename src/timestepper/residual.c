#include "timestepper.h"

PetscErrorCode computeResidual(SNES snes, 
                               Vec primPetscVec,
                               Vec residualPetscVec,
                               void *ptr)
{
  struct timeStepper *ts = (struct timeStepper*)ptr;

  int X1Start = ts->X1Start;
  int X2Start = ts->X2Start;
  int X1Size = ts->X1Size;
  int X2Size = ts->X2Size;

  if (ts->computeOldSourceTermsAndOldDivOfFluxes)
  {
    Vec primPetscVecOldLocal, primPetscVecHalfStepLocal;
    DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);
    DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecHalfStepLocal);

    /* Exchange ghost zone data. */
    if (ts->computeDivOfFluxAtTimeN)
    {
      /* Compute Div(flux) at t=n */
      DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                           ts->primPetscVecOld,
                           INSERT_VALUES,
                           primPetscVecOldLocal);
      DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                         ts->primPetscVecOld,
                         INSERT_VALUES,
                         primPetscVecOldLocal);
    }
    else if (ts->computeDivOfFluxAtTimeNPlusHalf)
    {
      /* Compute Div(flux) at t=n+1/2 */
      DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                           ts->primPetscVecHalfStep,
                           INSERT_VALUES,
                           primPetscVecHalfStepLocal);
      DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                         ts->primPetscVecHalfStep,
                         INSERT_VALUES,
                         primPetscVecHalfStepLocal);
    }

    REAL *primOldLocal, *primHalfStepLocal;
    REAL *divFluxOldGlobal, *sourceTermsOldGlobal, *conservedVarsOldGlobal;
    
    DMDAVecGetArray(ts->dmdaWithGhostZones,
                    primPetscVecOldLocal, &primOldLocal);

    DMDAVecGetArray(ts->dmdaWithGhostZones,
                    primPetscVecHalfStepLocal, &primHalfStepLocal);
    
    DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                    ts->divFluxPetscVecOld, &divFluxOldGlobal);

    DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                    ts->sourceTermsPetscVecOld, &sourceTermsOldGlobal); 
                    
    DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                    ts->conservedVarsPetscVecOld, &conservedVarsOldGlobal); 

    /* Loop through tiles. We use tiles to maximize cache usage.*/
    LOOP_OVER_TILES(X1Size, X2Size)
    {
      REAL primTile[TILE_SIZE];
      REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];

      /* Load data from the global memory on RAM onto a tile small enough
        * to reside on the cache */
      if (ts->computeDivOfFluxAtTimeN)
      {
        LOOP_INSIDE_TILE(-NG, TILE_SIZE_X2+NG, -NG, TILE_SIZE_X1+NG)
        {
          struct gridZone zone;
          setGridZone(iTile, jTile,
                      iInTile, jInTile,
                      X1Start, X2Start, 
                      X1Size, X2Size, 
                      &zone);
#pragma ivdep
          for (int var=0; var<DOF; var++)
          {
            primTile[INDEX_TILE(zone, var)] =
            primOldLocal[INDEX_LOCAL(zone, var)];
          }
        }
      } 
      else if (ts->computeDivOfFluxAtTimeNPlusHalf)
      {
        LOOP_INSIDE_TILE(-NG, TILE_SIZE_X2+NG, -NG, TILE_SIZE_X1+NG)
        {
          struct gridZone zone;
          setGridZone(iTile, jTile,
                      iInTile, jInTile,
                      X1Start, X2Start, 
                      X1Size, X2Size, 
                      &zone);
#pragma ivdep
          for (int var=0; var<DOF; var++)
          {
            primTile[INDEX_TILE(zone, var)] =
            primHalfStepLocal[INDEX_LOCAL(zone, var)];
          }
        }
      }

      /* Sync point */
      
      /* Apply boundary conditions on each tile */
      LOOP_INSIDE_TILE(0, TILE_SIZE_X2, 0, TILE_SIZE_X1)
      {
        struct gridZone zone;
        setGridZone(iTile, jTile,
                    iInTile, jInTile,
                    X1Start, X2Start, 
                    X1Size, X2Size, 
                    &zone);
        
        setZoneBoundaryFlags(&zone);
  
        applyTileBoundaryConditions(zone, primOldLocal, primTile);
      }

      /* Sync point */
  
      /* Work on the tiles.*/
      computeFluxesOverTile(primTile,
                            iTile, jTile,
                            X1Start, X2Start,
                            X1Size, X2Size,
                            fluxX1Tile, fluxX2Tile);

      LOOP_INSIDE_TILE(0, TILE_SIZE_X2, 0, TILE_SIZE_X1)
      {
        struct gridZone zone;
        setGridZone(iTile, jTile,
                    iInTile, jInTile,
                    X1Start, X2Start, 
                    X1Size, X2Size, 
                    &zone);

        REAL XCoords[NDIM], sourceTerms[DOF], conservedVars[DOF];

        getXCoords(&zone, CENTER, XCoords);
        struct geometry geom; setGeometry(XCoords, &geom);

        /* Now we need to compute conservedVarsOld using data from
         * primOldLocal. */
        struct fluidElement elem;
        setFluidElement(&primOldLocal[INDEX_LOCAL(zone, 0)], &geom, &elem);
        computeFluxes(&elem, &geom, 0, conservedVars);
        for (int var=0; var<DOF; var++)
        {
          conservedVarsOldGlobal[INDEX_GLOBAL(zone, var)] = 
            conservedVars[var];

          divFluxOldGlobal[INDEX_GLOBAL(zone, var)] = 
            (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(zone, var)]
             - fluxX1Tile[INDEX_TILE(zone, var)]
            )/zone.dX1
          + 
            (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(zone, var)]
             - fluxX2Tile[INDEX_TILE(zone, var)]
            )/zone.dX2;
        }

        if (ts->computeSourceTermsAtTimeN)
        {
          computeSourceTerms(&elem, &geom, XCoords, sourceTerms);

          for (int var=0; var<DOF; var++)
          {
            sourceTermsOldGlobal[INDEX_GLOBAL(zone, var)] =
              sourceTerms[var];
          }
        }
        else if (ts->computeSourceTermsAtTimeNPlusHalf)
        {
          setFluidElement(&primHalfStepLocal[INDEX_LOCAL(zone, 0)],
                          &geom, &elem);

          computeSourceTerms(&elem, &geom, XCoords, sourceTerms);

          for (int var=0; var<DOF; var++)
          {
            sourceTermsOldGlobal[INDEX_GLOBAL(zone, var)] =
              sourceTerms[var];
          }
        }

      }

    } /* End of LOOP_OVER_TILES */

    DMDAVecRestoreArray(ts->dmdaWithGhostZones,
                        primPetscVecOldLocal, &primOldLocal);

    DMDAVecRestoreArray(ts->dmdaWithGhostZones,
                        primPetscVecHalfStepLocal, &primHalfStepLocal);
    
    DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                        ts->divFluxPetscVecOld, &divFluxOldGlobal);

    DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                        ts->sourceTermsPetscVecOld, &sourceTermsOldGlobal); 
                    
    DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                        ts->conservedVarsPetscVecOld, &conservedVarsOldGlobal); 

    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);

    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecHalfStepLocal);

    /* All old sources and divFluxes have now been computed */
    ts->computeOldSourceTermsAndOldDivOfFluxes = 0;
  }

  /* The following computation requires no communication*/
  REAL *residualGlobal;
  REAL *divFluxOldGlobal, *sourceTermsOldGlobal, *conservedVarsOldGlobal;

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  residualPetscVec, &residualGlobal);

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  ts->divFluxPetscVecOld, &divFluxOldGlobal);

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  ts->sourceTermsPetscVecOld, &sourceTermsOldGlobal); 

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  ts->conservedVarsPetscVecOld, &conservedVarsOldGlobal); 

#if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
  REAL *primGlobal;

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  primPetscVec, &primGlobal);

#elif (TIME_STEPPING==IMPLICIT)
  Vec primPetscVecLocal;
  DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecLocal);

  /* Exchange ghost zone data. */
  DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                       ts->primPetscVec,
                       INSERT_VALUES,
                       primPetscVecLocal);
  DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                     ts->primPetscVec,
                     INSERT_VALUES,
                     primPetscVecLocal);

  REAL *primLocal;
  
  DMDAVecGetArray(ts->dmdaWithGhostZones,
                  primPetscVecOldLocal, &primOldLocal);
#endif

  LOOP_OVER_TILES(X1Size, X2Size)
  {
    REAL primTile[TILE_SIZE];

#if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
    LOOP_INSIDE_TILE(0, TILE_SIZE_X2, 0, TILE_SIZE_X1)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  X1Start, X2Start, 
                  X1Size, X2Size, 
                  &zone);
#pragma ivdep
      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(zone, var)] =
        primGlobal[INDEX_GLOBAL(zone, var)];
      }
    }
#elif (TIME_STEPPING==IMPLICIT)
    LOOP_INSIDE_TILE(-NG, TILE_SIZE_X2+NG, -NG, TILE_SIZE_X1+NG)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  X1Start, X2Start, 
                  X1Size, X2Size, 
                  &zone);
#pragma ivdep
      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(zone, var)] =
        primLocal[INDEX_LOCAL(zone, var)];
      }
    }
    /* Sync point */
    
    /* Apply boundary conditions on each tile */
    LOOP_INSIDE_TILE(0, TILE_SIZE_X2, 0, TILE_SIZE_X1)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  X1Start, X2Start, 
                  X1Size, X2Size, 
                  &zone);
      
      setZoneBoundaryFlags(&zone);
  
      applyTileBoundaryConditions(zone, primLocal, primTile);
    }

    REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];
    computeFluxesOverTile(primTile, 
                          iTile, jTile,
                          X1Start, X2Start,
                          X1Size, X2Size,
                          fluxX1Tile, fluxX2Tile);

#endif

    LOOP_INSIDE_TILE(0, TILE_SIZE_X2, 0, TILE_SIZE_X1)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  X1Start, X2Start, 
                  X1Size, X2Size, 
                  &zone);

      REAL XCoords[NDIM];

      getXCoords(&zone, CENTER, XCoords);
      struct geometry geom; setGeometry(XCoords, &geom);
      struct fluidElement elem;
      setFluidElement(&primTile[INDEX_TILE(zone, 0)], &geom, &elem);

      REAL conservedVars[DOF];
      computeFluxes(&elem, &geom, 0, conservedVars);

      for (int var=0; var<DOF; var++)
      {
#if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
        residualGlobal[INDEX_GLOBAL(zone, var)] = 
        (  conservedVars[var]
         - conservedVarsOldGlobal[INDEX_GLOBAL(zone, var)]
        )/ts->dt
      + divFluxOldGlobal[INDEX_GLOBAL(zone, var)]
      - sourceTermsOldGlobal[INDEX_GLOBAL(zone, var)]; 

#elif (TIME_STEPPING==IMPLICIT)
        residualGlobal[INDEX_GLOBAL(zone, var)] = 
        (  conservedVars[var]
         - conservedVarsOldGlobal[INDEX_GLOBAL(zone, var)]
        )/ts->dt
      + 0.5*(  divFluxOldGlobal[INDEX_GLOBAL(zone, var)]
             + 
               (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(zone, var)]
                - fluxX1Tile[INDEX_TILE(zone, var)]
               )/zone.dX1
              + 
               (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(zone, var)]
                - fluxX2Tile[INDEX_TILE(zone, var)]
               )/zone.dX2;
            )
      - sourceTermsOldGlobal[INDEX_GLOBAL(zone, var)]; 
#endif
      }
    }

  }

#if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
  DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                      primPetscVec, &primGlobal);

#elif (TIME_STEPPING==IMPLICIT)
  DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecLocal);

  DMDAVecRestoreArray(ts->dmdaWithGhostZones,
                      primPetscVecOldLocal, &primOldLocal);
#endif

  DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                      residualPetscVec, &residualGlobal);

  DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                      ts->divFluxPetscVecOld, &divFluxOldGlobal);

  DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                      ts->sourceTermsPetscVecOld, &sourceTermsOldGlobal); 

  DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                      ts->conservedVarsPetscVecOld, &conservedVarsOldGlobal); 

  return(0);
}

void computeFluxesOverTile(const REAL primTile[ARRAY_ARGS TILE_SIZE],
                           const int iTile, const int jTile,
                           const int X1Start, const int X2Start,
                           const int X1Size, const int X2Size,
                           REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
                           REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE])
{
  REAL fluxTileLeft[TILE_SIZE], fluxTileRight[TILE_SIZE];
  REAL conservedVarsTileLeft[TILE_SIZE], conservedVarsTileRight[TILE_SIZE];

  LOOP_INSIDE_TILE(-NG+1, TILE_SIZE_X2+NG-1, -NG+1, TILE_SIZE_X1+NG-1)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);
    
    REAL slope[DOF], primEdge[DOF];

    slopeLim(&primTile[INDEX_TILE_MINUS_ONE_X1(zone, 0)],
             &primTile[INDEX_TILE(zone, 0)],
             &primTile[INDEX_TILE_PLUS_ONE_X1(zone, 0)],
             slope);

    REAL XCoords[NDIM];
  
    /* Left Edge */
    for (int var=0; var<DOF; var++)
    {
      primEdge[var] = primTile[INDEX_TILE(zone, var)] - 0.5*slope[var];
    }
    getXCoords(&zone, FACE_X1, XCoords);
    struct geometry geom; setGeometry(XCoords, &geom);
    struct fluidElement elem;
    setFluidElement(primEdge, &geom, &elem);
  
    computeFluxes(&elem, &geom, 1,
                  &fluxTileLeft[INDEX_TILE(zone, 0)]);
    computeFluxes(&elem, &geom, 0,
                  &conservedVarsTileLeft[INDEX_TILE(zone, 0)]);

    /* Right Edge */
    for (int var=0; var<DOF; var++)
    {
      primEdge[var] = primTile[INDEX_TILE(zone, var)] + 0.5*slope[var];
    }
    getXCoords(&zone, FACE_X1_PLUS_ONE, XCoords);
    setGeometry(XCoords, &geom);
    setFluidElement(primEdge, &geom, &elem);
  
    computeFluxes(&elem, &geom, 1,
                  &fluxTileRight[INDEX_TILE(zone, 0)]);
    computeFluxes(&elem, &geom, 0,
                  &conservedVarsTileRight[INDEX_TILE(zone, 0)]);
  }
        
  LOOP_INSIDE_TILE(-NG+1, TILE_SIZE_X2+NG-1, -NG+1, TILE_SIZE_X1+NG-1)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

//    riemannSolver(fluxTileRight[INDEX_TILE_MINUS_ONE_X1(zone, 0)],
//                  fluxX1TileLeft[INDEX_TILE(zone, 0)],
//              conservedVarsTileRight[INDEX_TILE_MINUS_ONE_X1(zone, 0)],
//              conservedVarsTileLeft[INDEX_TILE(zone, 0)],
//              fluxX1Tile[INDEX_TILE(zone, 0)]);
  }

#if (COMPUTE_DIM==2)
  LOOP_INSIDE_TILE(-NG+1, TILE_SIZE_X2+NG-1, -NG+1, TILE_SIZE_X1+NG-1)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);
    
    REAL slope[DOF], primEdge[DOF];

    slopeLim(&primTile[INDEX_TILE_MINUS_ONE_X2(zone, 0)],
             &primTile[INDEX_TILE(zone, 0)],
             &primTile[INDEX_TILE_PLUS_ONE_X2(zone, 0)],
             slope);

    REAL XCoords[NDIM];
  
    /* Left Edge */
    for (int var=0; var<DOF; var++)
    {
      primEdge[var] = primTile[INDEX_TILE(zone, var)] - 0.5*slope[var];
    }
    getXCoords(&zone, FACE_X2, XCoords);
    struct geometry geom; setGeometry(XCoords, &geom);
    struct fluidElement elem;
    setFluidElement(primEdge, &geom, &elem);
  
    computeFluxes(&elem, &geom, 2,
                  &fluxTileLeft[INDEX_TILE(zone, 0)]);
    computeFluxes(&elem, &geom, 0,
                  &conservedVarsTileLeft[INDEX_TILE(zone, 0)]);

    /* Right Edge */
    for (int var=0; var<DOF; var++)
    {
      primEdge[var] = primTile[INDEX_TILE(zone, var)] + 0.5*slope[var];
    }
    getXCoords(&zone, FACE_X2_PLUS_ONE, XCoords);
    setGeometry(XCoords, &geom);
    setFluidElement(primEdge, &geom, &elem);
  
    computeFluxes(&elem, &geom, 2,
                  &fluxTileRight[INDEX_TILE(zone, 0)]);
    computeFluxes(&elem, &geom, 0,
                  &conservedVarsTileRight[INDEX_TILE(zone, 0)]);

  }
        
  LOOP_INSIDE_TILE(-NG+1, TILE_SIZE_X2+NG-1, -NG+1, TILE_SIZE_X1+NG-1)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

//    riemannSolver(fluxTileRight[INDEX_TILE_MINUS_ONE_X2(zone, 0)],
//                  fluxX1TileLeft[INDEX_TILE(zone, 0)],
//              conservedVarsTileRight[INDEX_TILE_MINUS_ONE_X2(zone, 0)],
//              conservedVarsTileLeft[INDEX_TILE(zone, 0)],
//              fluxX2Tile[INDEX_TILE(zone, 0)]);
  }
#endif /* COMPUTE_DIM==2 */

}
