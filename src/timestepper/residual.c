#include "timestepper.h"

PetscErrorCode computeResidual(SNES snes, 
                               Vec primPetscVec,
                               Vec residualPetscVec,
                               void *ptr)
{
  struct timeStepper *ts = (struct timeStepper*)ptr;

  int X1Start, X2Start, X3Start;
  int X1Size, X2Size, X3Size;

  DMDAGetCorners(ts->dmdaWithGhostZones,
                 &X1Start, &X2Start, &X3Start,
                 &X1Size, &X2Size, &X3Size);

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

    ARRAY(primOldLocal);
    ARRAY(primHalfStepLocal);
    ARRAY(divFluxOldGlobal);
    ARRAY(sourceTermsOldGlobal);
    ARRAY(conservedVarsOldGlobal);
    
    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecOldLocal,
                       &primOldLocal);
    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecHalfStepLocal,
                       &primHalfStepLocal);
    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->divFluxPetscVecOld,
                       &divFluxOldGlobal);
    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->sourceTermsPetscVecOld,
                       &sourceTermsOldGlobal);           
    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->conservedVarsPetscVecOld,
                       &conservedVarsOldGlobal); 

    /* Loop through tiles. We use tiles to maximize cache usage.*/
    LOOP_OVER_TILES(X1Size, X2Size)
    {
      REAL primTile[TILE_SIZE];
      REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];

      /* Load data from the global memory on RAM onto a tile small enough
        * to reside on the cache */
      if (ts->computeDivOfFluxAtTimeN)
      {
        LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
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
            primTile[INDEX_TILE(&zone, var)] =
            INDEX_PETSC(primOldLocal, &zone, var);
          }
        }
      } 
      else if (ts->computeDivOfFluxAtTimeNPlusHalf)
      {
        LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
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
            primTile[INDEX_TILE(&zone, var)] =
            INDEX_PETSC(primHalfStepLocal, &zone, var);
          }
        }
      }

      /* Sync point */
      
//      /* Apply boundary conditions on each tile */
//      LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
//      {
//        struct gridZone zone;
//        setGridZone(iTile, jTile,
//                    iInTile, jInTile,
//                    X1Start, X2Start, 
//                    X1Size, X2Size, 
//                    &zone);
//        
//        setZoneBoundaryFlags(&zone);
//  
//        applyTileBoundaryConditions(&zone, primOldLocal, primTile);
//      }

      /* Sync point */
  
      /* Work on the tiles.*/
      computeFluxesOverTile(primTile,
                            iTile, jTile,
                            X1Start, X2Start,
                            X1Size, X2Size,
                            fluxX1Tile, fluxX2Tile);

      LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
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
        setFluidElement(&INDEX_PETSC(primOldLocal, &zone, 0), &geom, &elem);
        computeFluxes(&elem, &geom, 0, conservedVars);
        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(conservedVarsOldGlobal, &zone, var) = 
            conservedVars[var];

          INDEX_PETSC(divFluxOldGlobal, &zone, var) = 
            (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(&zone, var)]
             - fluxX1Tile[INDEX_TILE(&zone, var)]
            )/zone.dX1
          + 
            (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(&zone, var)]
             - fluxX2Tile[INDEX_TILE(&zone, var)]
            )/zone.dX2;
        }

        if (ts->computeSourceTermsAtTimeN)
        {
          computeSourceTerms(&elem, &geom, XCoords, sourceTerms);

          for (int var=0; var<DOF; var++)
          {
            INDEX_PETSC(sourceTermsOldGlobal, &zone, var) = 
              sourceTerms[var];
          }
        }
        else if (ts->computeSourceTermsAtTimeNPlusHalf)
        {
          setFluidElement(&INDEX_PETSC(primHalfStepLocal, &zone, 0),
                          &geom, &elem);

          computeSourceTerms(&elem, &geom, XCoords, sourceTerms);

          for (int var=0; var<DOF; var++)
          {
            INDEX_PETSC(sourceTermsOldGlobal, &zone, var) = 
              sourceTerms[var];
          }
        }

      }

    } /* End of LOOP_OVER_TILES */

    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecOldLocal,
                           &primOldLocal);
    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecHalfStepLocal,
                           &primHalfStepLocal);
    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->divFluxPetscVecOld,
                           &divFluxOldGlobal);
    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->sourceTermsPetscVecOld,
                           &sourceTermsOldGlobal);           
    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->conservedVarsPetscVecOld,
                           &conservedVarsOldGlobal); 

    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);
    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecHalfStepLocal);

    /* All old sources and divFluxes have now been computed */
    ts->computeOldSourceTermsAndOldDivOfFluxes = 0;
  }

  /* The following computation requires no communication*/
  ARRAY(residualGlobal);
  ARRAY(divFluxOldGlobal);
  ARRAY(sourceTermsOldGlobal);
  ARRAY(conservedVarsOldGlobal);

  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, residualPetscVec,
                     &residualGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->divFluxPetscVecOld,
                     &divFluxOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->sourceTermsPetscVecOld,
                     &sourceTermsOldGlobal); 
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->conservedVarsPetscVecOld,
                     &conservedVarsOldGlobal); 

  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)

    ARRAY(primGlobal);
    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, primPetscVec, &primGlobal);

  #elif (TIME_STEPPING==IMPLICIT)
    Vec primPetscVecLocal;
    DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecLocal);

    /* Exchange ghost zone data. */
    DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                         primPetscVec,
                         INSERT_VALUES,
                         primPetscVecLocal);
    DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                       primPetscVec,
                       INSERT_VALUES,
                       primPetscVecLocal);

    ARRAY(primLocal);
    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecLocal, &primLocal);
  #endif

  LOOP_OVER_TILES(X1Size, X2Size)
  {
    #if (TIME_STEPPING==IMPLICIT)
      REAL primTile[TILE_SIZE];

      LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
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
          primTile[INDEX_TILE(&zone, var)] =
          INDEX_PETSC(primLocal, &zone, var);
        }
      }
      /* Sync point */
    
      /* Apply boundary conditions on each tile */
//      LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
//      {
//        struct gridZone zone;
//        setGridZone(iTile, jTile,
//                    iInTile, jInTile,
//                    X1Start, X2Start, 
//                    X1Size, X2Size, 
//                    &zone);
//      
//        setZoneBoundaryFlags(&zone);
//  
//        applyTileBoundaryConditions(&zone, primLocal, primTile);
//      }

      REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];
      computeFluxesOverTile(primTile, 
                            iTile, jTile,
                            X1Start, X2Start,
                            X1Size, X2Size,
                            fluxX1Tile, fluxX2Tile);

    #endif

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
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
      #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
        setFluidElement(&INDEX_PETSC(primGlobal, &zone, 0), &geom, &elem);
      #elif (TIME_STEPPING==IMPLICIT)
        setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
      #endif

      REAL conservedVars[DOF];
      computeFluxes(&elem, &geom, 0, conservedVars);

      #if (TIME_STEPPING==IMEX)
        REAL sourceTerms[DOF];
        computeSourceTerms(&elem, &geom, XCoords, sourceTerms);
      #endif

      for (int var=0; var<DOF; var++)
      {
        #if (TIME_STEPPING==EXPLICIT)

          INDEX_PETSC(residualGlobal, &zone, var) = 
          (  conservedVars[var]
           - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
          )/ts->dt
          + INDEX_PETSC(divFluxOldGlobal, &zone, var)
          - INDEX_PETSC(sourceTermsOldGlobal, &zone, var); 

        #elif (TIME_STEPPING==IMEX)

          INDEX_PETSC(residualGlobal, &zone, var) = 
          (  conservedVars[var]
           - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
          )/ts->dt
          + INDEX_PETSC(divFluxOldGlobal, &zone, var)
          - 0.5*(  INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
                 + sourceTerms[var]
                ); 

        #elif (TIME_STEPPING==IMPLICIT)

          INDEX_PETSC(residualGlobal, &zone, var) = 
          (  conservedVars[var]
          - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
          )/ts->dt
          + 0.5*(  INDEX_PETSC(divFluxOldGlobal, &zone, var)
                 + 
                  (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(&zone, var)]
                   - fluxX1Tile[INDEX_TILE(&zone, var)]
                  )/zone.dX1
                 + 
                  (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(&zone, var)]
                   - fluxX2Tile[INDEX_TILE(&zone, var)]
                  )/zone.dX2
                )
          - INDEX_PETSC(sourceTermsOldGlobal, &zone, var); 

        #endif
      }
    }

  }

  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, primPetscVec,
                           &primGlobal);

  #elif (TIME_STEPPING==IMPLICIT)
    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecLocal);

    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecLocal,
                           &primLocal);
  #endif

  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, residualPetscVec,
                         &residualGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->divFluxPetscVecOld,
                         &divFluxOldGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->sourceTermsPetscVecOld,
                         &sourceTermsOldGlobal); 
  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->conservedVarsPetscVecOld,
                         &conservedVarsOldGlobal); 
  return(0);
}
