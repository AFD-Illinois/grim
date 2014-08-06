#include "timestepper.h"

PetscErrorCode computeResidualIMEX(SNES snes, 
                                   Vec primPetscVec,
                                   Vec residualPetscVec,
                                   void *ptr)
{
  struct timeStepper *ts = (struct timeStepper*)ptr;

  int x1Start, x1Size;
  int x2Start, x2Size;
  int x3Start, x3Size;

  /* dmdaWithGhostZones and dmdaWithoutGhostZones will both give the same
   * start indices and sizes */
  DMDAGetCorners(ts->dmdaWithGhostZones,
                 &x1Start, &x2Start, &x3Start,
                 &x1Size, &x2Size, &x3Size);

  if (!ts->computedOldSourceTermsAndOldFluxes)
  {
    Vec primPetscVecOldLocal;
    DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);

    /* Exchange ghost zone data. */
    DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld,
                         INSERT_VALUES,
                         primPetscVecOldLocal);
    DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                       ts->primPetscVecOld,
                       INSERT_VALUES,
                       primPetscVecOldLocal);

    REAL *primOldLocal;
    REAL *divFluxOldGlobal, *sourceTermsOldGlobal, *conservedVarsOldGlobal;
    
    DMDAVecGetArray(ts->dmdaWithGhostZones,
                    primPetscVecOldLocal, &primOldLocal);
    
    DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                    ts->divFluxPetscVecOld, &divFluxOldGlobal);

    DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                    ts->sourceTermsPetscVecOld, &sourceTermsOldGlobal); 
                    
    DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                    ts->conservedVarsPetscVecOld, &conservedVarsOldGlobal); 

    struct gridZone zone;
    struct geometry geom;
    struct fluidElement elem;

    /* Loop through tiles. We use tiles to maximize cache usage.*/
#if (COMPUTE_DIM==2)
    for (int jTile=0; jTile<x2Size/TILE_SIZE_X2; jTile++)
    {
#else
    {
      int jTile=0; int jInTile=0; int x2Start=0; int x2Size = 0;
#endif
      for (int iTile=0; iTile<x1Size/TILE_SIZE_X1; iTile++)
      {
        REAL primTile[TILE_SIZE];
        REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];
        REAL fluxTileLeft[TILE_SIZE], fluxTileRight[TILE_SIZE];
        REAL conservedVarsTileLeft[TILE_SIZE], conservedVarsTileRight[TILE_SIZE];

      /* Load data from the global memory on RAM onto a tile small enough
       * to reside on the cache */
#if (COMPUTE_DIM==2)
        for (int jInTile=-NG; jInTile<TILE_SIZE_X2+NG; jInTile++)
#endif
          for (int iInTile=-NG; iInTile<TILE_SIZE_X1+NG; iInTile++)
          {
            setGridZone(iTile, jTile,
                        iInTile, jInTile,
                        x1Start, x2Start, 
                        x1Size, x2Size, 
                        &zone);
#pragma ivdep
            for (int var=0; var<DOF; var++)
            {
              primTile[INDEX_TILE(zone, var)] =
              primOldLocal[INDEX_LOCAL(zone, var)];
            }
          }

        /* Sync point */
      
        /* Apply boundary conditions on each tile */
#if (COMPUTE_DIM==2)
        for (int jInTile=0; jInTile<TILE_SIZE_X2; jInTile++)
#endif
          for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
          {
            setGridZone(iTile, jTile,
                        iInTile, jInTile,
                        x1Start, x2Start, 
                        x1Size, x2Size, 
                        &zone);
        
            setZoneBoundaryFlags(&zone);
  
            applyTileBoundaryConditions(&zone, primOldLocal, primTile);
          }

        /* Sync point */
  
        /* Work on the tiles.*/
#if (COMPUTE_DIM==2)
        for (int jInTile=-NG+1; jInTile<TILE_SIZE_X2+NG-1; jInTile++)
#endif
          for (int iInTile=-NG+1; iInTile<TILE_SIZE_X1+NG-1; iInTile++)
          {
            setGridZone(iTile, jTile,
                        iInTile, jInTile,
                        x1Start, x2Start, 
                        x1Size, x2Size, 
                        &zone);
            
            REAL slope[DOF], primEdge[DOF];

            slopeLim(primTile[INDEX_TILE_MINUS_ONE_X1(zone, 0)],
                     primTile[INDEX_TILE(zone, 0)],
                     primTile[INDEX_TILE_PLUS_ONE_X1(zone, 0)],
                     slope);

            REAL XCoords[NDIM];
  
            /* Left Edge */
            for (int var=0; var<DOF; var++)
            {
              primEdge[var] = primTile[INDEX_TILE(zone, var)] - 0.5*slope[var];
            }
            getXCoords(zone, FACE_X1, XCoords);
            setGeometry(XCoords, &geom);
            setFluidElement(primEdge, &geom, &elem);
  
            computeFluxesAndConservedVars(&elem, &geom, 1,
                          &fluxTileLeft[INDEX_TILE(zone, 0)],
                          &conservedVarsTileLeft[INDEX_TILE(zone, 0)]);

            /* Right Edge */
            for (int var=0; var<DOF; var++)
            {
              primEdge[var] = primTile[INDEX_TILE(zone, var)] + 0.5*slope[var];
            }
            getXCoords(zone, FACE_X1_PLUS_ONE, XCoords);
            setGeometry(XCoords, &geom);
            setFluidElement(primEdge, &geom, &elem);
  
            computeFluxesAndConservedVars(&elem, &geom, 1,
                          &fluxTileRight[INDEX_TILE(zone, 0)],
                          &conservedVarsTileRight[INDEX_TILE(zone, 0)]);

          }

        
#if (COMPUTE_DIM==2)
        for (int jInTile=-NG+1; jInTile<TILE_SIZE_X2+NG-1; jInTile++)
#endif
          for (int iInTile=-NG+1; iInTile<TILE_SIZE_X1+NG-1; iInTile++)
          {
            setGridZone(iTile, jTile,
                        iInTile, jInTile,
                        x1Start, x2Start, 
                        x1Size, x2Size, 
                        &zone);

            riemannSolver(fluxTileRight[INDEX_TILE_MINUS_ONE_X1(zone, 0)],
                          fluxX1TileLeft[INDEX_TILE(zone, 0)],
                      conservedVarsTileRight[INDEX_TILE_MINUS_ONE_X1(zone, 0)],
                      conservedVarsTileLeft[INDEX_TILE(zone, 0)],
                      fluxX1Tile[INDEX_TILE(zone, var)]);
          }

#if (COMPUTE_DIM==2)
        for (int jInTile=-NG+1; jInTile<TILE_SIZE_X2+NG-1; jInTile++)
        {
          for (int iInTile=-NG+1; iInTile<TILE_SIZE_X1+NG-1; iInTile++)
          {
            setGridZone(iTile, jTile,
                        iInTile, jInTile,
                        x1Start, x2Start, 
                        x1Size, x2Size, 
                        &zone);
            
            REAL slope[DOF], primEdge[DOF];

            slopeLim(primTile[INDEX_TILE_MINUS_ONE_X2(zone, 0)],
                     primTile[INDEX_TILE(zone, 0)],
                     primTile[INDEX_TILE_PLUS_ONE_X2(zone, 0)],
                     slope);

            REAL XCoords[NDIM];
  
            /* Left Edge */
            for (int var=0; var<DOF; var++)
            {
              primEdge[var] = primTile[INDEX_TILE(zone, var)] - 0.5*slope[var];
            }
            getXCoords(zone, FACE_X2, XCoords);
            setGeometry(XCoords, &geom);
            setFluidElement(primEdge, &geom, &elem);
  
            computeFluxesAndConservedVars(&elem, &geom, 2,
                          &fluxTileLeft[INDEX_TILE(zone, 0)],
                          &conservedVarsTileLeft[INDEX_TILE(zone, 0)]);

            /* Right Edge */
            for (int var=0; var<DOF; var++)
            {
              primEdge[var] = primTile[INDEX_TILE(zone, var)] + 0.5*slope[var];
            }
            getXCoords(zone, FACE_X2_PLUS_ONE, XCoords);
            setGeometry(XCoords, &geom);
            setFluidElement(primEdge, &geom, &elem);
  
            computeFluxesAndConservedVars(&elem, &geom, 2,
                          &fluxTileRight[INDEX_TILE(zone, 0)],
                          &conservedVarsTileRight[INDEX_TILE(zone, 0)]);

          }
        }
        
        for (int jInTile=-NG+1; jInTile<TILE_SIZE_X2+NG-1; jInTile++)
        {
          for (int iInTile=-NG+1; iInTile<TILE_SIZE_X1+NG-1; iInTile++)
          {
            setGridZone(iTile, jTile,
                        iInTile, jInTile,
                        x1Start, x2Start, 
                        x1Size, x2Size, 
                        &zone);

            riemannSolver(fluxTileRight[INDEX_TILE_MINUS_ONE_X2(zone, 0)],
                          fluxX1TileLeft[INDEX_TILE(zone, 0)],
                      conservedVarsTileRight[INDEX_TILE_MINUS_ONE_X2(zone, 0)],
                      conservedVarsTileLeft[INDEX_TILE(zone, 0)],
                      fluxX2Tile[INDEX_TILE(zone, var)]);
          }
        }
#endif /* COMPUTE_DIM==2 */

#if (COMPUTE_DIM==2)
        for (int jInTile=0; jInTile<TILE_SIZE_X2; jInTile++)
#endif
          for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
          {
            setGridZone(iTile, jTile,
                        iInTile, jInTile,
                        x1Start, x2Start, 
                        x1Size, x2Size, 
                        &zone);

            REAL XCoords[NDIM], sourceTerms[DOF];

            getXCoords(zone, CENTER, XCoords);
            setGeometry(XCoords, &geom);
            setFluidElement(primTile[INDEX_TILE(zone, 0)], &geom, &elem);

            /* Using sourceTerms in place of the fluxes argument cause we do not
             * need fluxes and I don't want to create fluxes[DOF] */
            computeFluxesAndConservedVars(&elem, &geom, 0,
                                          sourceTerms, conservedVars);

            computeSourceTerms(&elem, &geom, XCoords, sourceTerms);

            for (int var=0; var<DOF; var++)
            {
              divFluxOldGlobal[INDEX_GLOBAL(zone, var)] = 
                (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(zone, var)]
                 - fluxX1Tile[INDEX_TILE(zone, var)]
                )/zone->dX1
              + 
                (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(zone, var)]
                 - fluxX2Tile[INDEX_TILE(zone, var)]
                )/zone->dX2;

              sourceTermsOldGlobal[INDEX_GLOBAL(zone, var)] = 
                sourceTerms[var];

              conservedVarsOldGlobal[INDEX_GLOBAL(zone, var)] = 
                conservedVars[var];
            }
          }

      } /* End of iTile loop */
    } /* End of jTile loop */

    DMDAVecRestoreArray(ts->dmdaWithGhostZones,
                        primPetscVecOldLocal, &primOldLocal);
    
    DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                        ts->divFluxPetscVecOld, &divFluxOldGlobal);

    DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                        ts->sourceTermsPetscVecOld, &sourceTermsOldGlobal); 
                    
    DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                        ts->conservedVarsPetscVecOld, &conservedVarsOldGlobal); 

    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);

    /* All old sources and divFluxes have now been computed */
    ts->computedOldSourceTermsAndOldFluxes = 1;
  }

  /* The following computation requires no communication*/
  REAL *primGlobal, *residualGlobal;
  REAL *divFluxOldGlobal, *sourceTermsOldGlobal, *conservedVarsOldGlobal;

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  primPetscVec, &primGlobal);

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  residualPetscVec, &residualGlobal);

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  ts->divFluxPetscVecOld, &divFluxOldGlobal);

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  ts->sourceTermsPetscVecOld, &sourceTermsOldGlobal); 

  DMDAVecGetArray(ts->dmdaWithoutGhostZones,
                  ts->conservedVarsPetscVecOld, &conservedVarsOldGlobal); 

  struct gridZone zone;
  struct geometry geom;
  struct fluidElement elem;

  for (int iTile=0; iTile<x1Size/TILE_SIZE_X1; iTile++)
  {
    REAL primTile[TILE_SIZE];
    REAL conservedVarsTile[TILE_SIZE];

    for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
    {
      setGridZone1D(iTile, iInTile, x1Start, &zone);
#pragma ivdep
      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(zone, var)] =
        primGlobal[INDEX_GLOBAL(zone, var)];
      }
    }

    for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
    {
      setGridZone1D(iTile, iInTile, x1Start, &zone);

      REAL XCoords[NDIM];

      getXCoords1D(zone, CENTER, XCoords);
      setGeometry(XCoords, &geom);
      setFluidElement(&primTile[INDEX_TILE(zone, 0)], &geom, &elem);

      computeFluxes(&elem, &geom, 0,
                    &conservedVarsTile[INDEX_TILE(zone, 0)]);
    }

    for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
    {
      setGridZone1D(iTile, iInTile, x1Start, &zone);

      for (int var=0; var<DOF; var++)
      {
#if (PROBLEM==CONS_TO_PRIM_INVERSION_TEST)
        /* Performing the CONS_TO_PRIM_INVERSION_TEST */
        residualGlobal[INDEX_GLOBAL(zone, var)] = 
        conservedVarsTile[INDEX_TILE(zone, var)] - 
        conservedVarsOldGlobal[INDEX_GLOBAL(zone, var)];
#else
        /* Production mode */
        residualGlobal[INDEX_GLOBAL(zone, var)] = 
        (conservedVarsTile[INDEX_TILE(zone, var)] -
         conservedVarsOldGlobal[INDEX_GLOBAL(zone, var)])/ts->dt
      + divFluxOldGlobal[INDEX_GLOBAL(zone, var)]
      + sourceTermsOldGlobal[INDEX_GLOBAL(zone, var)]; 
#endif
      }
    }

  }

  DMDAVecRestoreArray(ts->dmdaWithoutGhostZones,
                      primPetscVec, &primGlobal);

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
