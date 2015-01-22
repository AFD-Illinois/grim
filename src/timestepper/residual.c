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

  ARRAY(primGlobal);
  ARRAY(primOldGlobal);
  ARRAY(primHalfStepGlobal);
  ARRAY(divFluxOldGlobal);
  ARRAY(sourceTermsOldGlobal);
  ARRAY(conservedVarsOldGlobal);
  ARRAY(connectionGlobal);
  ARRAY(dtGlobal);
  ARRAY(residualGlobal);

  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, primPetscVec, 
                     &primGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecOld,
                     &primOldGlobal); 
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
                     &primHalfStepGlobal); 

  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->divFluxPetscVecOld,
                     &divFluxOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->sourceTermsPetscVecOld,
                     &sourceTermsOldGlobal);           
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->conservedVarsPetscVecOld,
                     &conservedVarsOldGlobal); 
  DMDAVecGetArrayDOF(ts->connectionDMDA, ts->connectionPetscVec,
                     &connectionGlobal);

  DMDAVecGetArrayDOF(ts->dmdaDt, ts->dtPetscVec, &dtGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, residualPetscVec,
                     &residualGlobal);

  #if (CONDUCTION)
    ARRAY(gradTGlobal);
    ARRAY(graduConGlobal);
    ARRAY(graduConHigherOrderTerm1Global);

    DMDAVecGetArrayDOF(ts->gradTDM, ts->gradTPetscVec, &gradTGlobal);
    DMDAVecGetArrayDOF(ts->graduConDM, ts->graduConPetscVec, 
                       &graduConGlobal);
    DMDAVecGetArrayDOF(ts->graduConHigherOrderTerm1DM, 
                       ts->graduConHigherOrderTerm1PetscVec,
                       &graduConHigherOrderTerm1Global);
  #endif

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

    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecOldLocal,
                       &primOldLocal);
    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecHalfStepLocal,
                       &primHalfStepLocal);


    /* Loop through tiles. We use tiles to maximize cache usage.*/
    #if (USE_OPENMP)
      #pragma omp parallel for
    #endif
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
          for (int var=0; var<DOF; var++)
          {
            primTile[INDEX_TILE(&zone, var)] =
            INDEX_PETSC(primHalfStepLocal, &zone, var);
          }
        }
      }

      /* Sync point */
      
      applyTileBoundaryConditions(iTile, jTile,
                                  X1Start, X2Start,
                                  X1Size, X2Size,
                                  primTile);

      applyAdditionalProblemSpecificBCs(iTile, jTile,
                                        X1Start, X2Start,
                                        X1Size, X2Size,
                                        primTile);

      /* Sync point */
  
      /* Work on the tiles.*/
      computeFluxesOverTile(primTile,
                            iTile, jTile,
                            X1Start, X2Start,
                            X1Size, X2Size,
                            fluxX1Tile, fluxX2Tile,
                            dtGlobal);

      applyProblemSpecificFluxFilter(iTile, jTile,
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
          #if (COMPUTE_DIM==2)
          + 
            (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(&zone, var)]
             - fluxX2Tile[INDEX_TILE(&zone, var)]
            )/zone.dX2
          #endif
            ;
        }

        if (ts->computeSourceTermsAtTimeN)
        {
          computeSourceTerms(&elem, &geom,
                             &INDEX_PETSC(connectionGlobal, &zone, 0),
                             sourceTerms);

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

          computeSourceTerms(&elem, &geom,
                             &INDEX_PETSC(connectionGlobal, &zone, 0),
                             sourceTerms);

          for (int var=0; var<DOF; var++)
          {
            INDEX_PETSC(sourceTermsOldGlobal, &zone, var) = 
              sourceTerms[var];
          }
        }

      }

      #if (CONDUCTION)
        addConductionSourceTermsToResidual
          (primTile,
           primGlobal, primHalfStepGlobal, primOldGlobal,
           connectionGlobal, 
           gradTGlobal, graduConGlobal, graduConHigherOrderTerm1Global,
           ts->dt,
           ts->computeOldSourceTermsAndOldDivOfFluxes,
           ts->computeDivOfFluxAtTimeN,
           ts->computeDivOfFluxAtTimeNPlusHalf,
           iTile, jTile, X1Start, X2Start, X1Size, X2Size,
           residualGlobal
          );
      #endif

    } /* End of LOOP_OVER_TILES */

    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecOldLocal,
                           &primOldLocal);
    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecHalfStepLocal,
                           &primHalfStepLocal);

    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);
    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecHalfStepLocal);

    /* All old sources and divFluxes have now been computed */
    ts->computeOldSourceTermsAndOldDivOfFluxes = 0;
  }

  /* The following computation requires no communication*/

  #if (TIME_STEPPING==IMPLICIT)
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

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
  LOOP_OVER_TILES(X1Size, X2Size)
  {
    REAL primTile[TILE_SIZE];
    REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];

    #if (TIME_STEPPING==IMPLICIT)
      LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
      {
        struct gridZone zone;
        setGridZone(iTile, jTile,
                    iInTile, jInTile,
                    X1Start, X2Start, 
                    X1Size, X2Size, 
                    &zone);
        for (int var=0; var<DOF; var++)
        {
          primTile[INDEX_TILE(&zone, var)] =
          INDEX_PETSC(primLocal, &zone, var);
        }
      }
      /* Sync point */
    
      /* Apply boundary conditions on each tile */
      applyTileBoundaryConditions(iTile, jTile,
                                  X1Start, X2Start,
                                  X1Size, X2Size,
                                  primTile);

      applyAdditionalProblemSpecificBCs(iTile, jTile,
                                        X1Start, X2Start,
                                        X1Size, X2Size,
                                        primTile);

      computeFluxesOverTile(primTile, 
                            iTile, jTile,
                            X1Start, X2Start,
                            X1Size, X2Size,
                            fluxX1Tile, fluxX2Tile,
                            dtGlobal);

      applyProblemSpecificFluxFilter(iTile, jTile,
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
      setFluidElement(&INDEX_PETSC(primGlobal, &zone, 0), &geom, &elem);

      REAL conservedVars[DOF];
      computeFluxes(&elem, &geom, 0, conservedVars);

      #if (TIME_STEPPING==IMEX || TIME_STEPPING==IMPLICIT)
        REAL sourceTerms[DOF];
        computeSourceTerms(&elem, &geom,
                           &INDEX_PETSC(connectionGlobal, &zone, 0),
                           sourceTerms);
      #endif

      REAL g = sqrt(-geom.gDet);
      REAL norm = g;

      for (int var=0; var<DOF; var++)
      {
        #if (TIME_STEPPING==EXPLICIT)

          INDEX_PETSC(residualGlobal, &zone, var) = 
          ( (  conservedVars[var]
             - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
            )/ts->dt
            + INDEX_PETSC(divFluxOldGlobal, &zone, var)
            - INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
          )/norm;

        #elif (TIME_STEPPING==IMEX)

          INDEX_PETSC(residualGlobal, &zone, var) = 
          ( (  conservedVars[var]
             - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
            )/ts->dt
            + INDEX_PETSC(divFluxOldGlobal, &zone, var)
            - 0.5*(  INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
                   + sourceTerms[var]
                  )
          )/norm;

        #elif (TIME_STEPPING==IMPLICIT)

          INDEX_PETSC(residualGlobal, &zone, var) = 
          ( (  conservedVars[var]
             - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
            )/ts->dt
            + 0.5*(  INDEX_PETSC(divFluxOldGlobal, &zone, var)
                   + 
                    (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(&zone, var)]
                     - fluxX1Tile[INDEX_TILE(&zone, var)]
                    )/zone.dX1
                  #if (COMPUTE_DIM==2)
                   + 
                    (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(&zone, var)]
                     - fluxX2Tile[INDEX_TILE(&zone, var)]
                    )/zone.dX2
                  #endif
                  )
            - 0.5*(  INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
                   + sourceTerms[var]
                  )
          )/norm;

        #endif
      }

    }

    #if (CONDUCTION)
      addConductionSourceTermsToResidual
        (primTile,
         primGlobal, primHalfStepGlobal, primOldGlobal,
         connectionGlobal,
         gradTGlobal, graduConGlobal, graduConHigherOrderTerm1Global,
         ts->dt,
         ts->computeOldSourceTermsAndOldDivOfFluxes,
         ts->computeDivOfFluxAtTimeN,
         ts->computeDivOfFluxAtTimeNPlusHalf,
         iTile, jTile, X1Start, X2Start, X1Size, X2Size,
         residualGlobal
        );
    #endif


  } /* End of LOOP_OVER_TILES */


  #if (TIME_STEPPING==IMPLICIT)
    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecLocal);

    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecLocal,
                           &primLocal);
  #endif

  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, primPetscVec,
                         &primGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecOld,
                         &primOldGlobal); 
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
                         &primHalfStepGlobal);

  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->divFluxPetscVecOld,
                         &divFluxOldGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->sourceTermsPetscVecOld,
                         &sourceTermsOldGlobal); 
  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->conservedVarsPetscVecOld,
                         &conservedVarsOldGlobal); 
  DMDAVecRestoreArrayDOF(ts->connectionDMDA, ts->connectionPetscVec,
                         &connectionGlobal);

  DMDAVecRestoreArrayDOF(ts->dmdaDt, ts->dtPetscVec, &dtGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, residualPetscVec,
                         &residualGlobal);
  #if (CONDUCTION)
    DMDAVecRestoreArrayDOF(ts->gradTDM, ts->gradTPetscVec, &gradTGlobal);
    DMDAVecRestoreArrayDOF(ts->graduConDM, ts->graduConPetscVec, 
                           &graduConGlobal);
    DMDAVecRestoreArrayDOF(ts->graduConHigherOrderTerm1DM, 
                           ts->graduConHigherOrderTerm1PetscVec,
                           &graduConHigherOrderTerm1Global);
  #endif
  return(0);
}
