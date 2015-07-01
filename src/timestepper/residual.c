#include "timestepper.h"

PetscErrorCode computeResidual(SNES snes, 
                               Vec primVec,
                               Vec residualVec,
                               void *ptr
                              )
{
  struct timeStepper *ts = (struct timeStepper*)ptr;
  struct gridData *prim  = &ts->primN;

  setPointerToVec(&ts->connection);
  setPointerToVec(&ts->dtGrid);

  if (ts->isZerothIterationOfSNES)
  {
    if (ts->computeDivOfFluxAtNPlusHalf)
    {
      /* Replace prim = &ts->primN, with prim = &ts->primNPlusHalf */
      prim = &ts->primNPlusHalf;
    }

    startFillingVecGhost(prim);

    setPointerToVec(&ts->primN);
    setPointerToVec(&ts->primNPlusHalf);

    if (ts->computeConsVarsAndSourcesAtN)
    {
      LOOP_OVER_TILES(&ts->primN)
      {
        struct gridTile tile;
        setTile(iTile, jTile, kGlobal, &ts->primN, &tile);

        LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
        {
          struct gridZone zone;
          setZone(iInTile, jInTile, 0, &tile, &zone);

          REAL XCoords[NDIM];
          REAL primVars[DOF], sources[DOF], conservedVars[DOF];
          REAL connectionCoeffs[64];

          getXCoords(&zone, CENTER, XCoords);
          struct geometry geom; setGeometry(XCoords, &geom);

          memcpy(primVars, &INDEX_GRID(&ts->primN, &zone, 0),
                 sizeof(REAL[DOF])
                );

          struct fluidElement elem;
          setFluidElement(primVars, &geom, &elem);
          computeFluxes(&elem, &geom, 0, conservedVars);

          memcpy(connectionCoeffs, &INDEX_GRID(&ts->connection, &zone, 0),
                 sizeof(REAL[64])
                );

          computeSourceTerms(&elem, &geom, connectionCoeffs, sources);

          memcpy(&INDEX_GRID(&ts->conservedVarsN, &zone, 0), conservedVars,
                 sizeof(REAL[DOF])
                );

          memcpy(&INDEX_GRID(&ts->sourcesN, &zone, 0), sources,
                 sizeof(REAL[DOF])
                );
        }
      }
    }

    if (ts->computeSourcesAtNPlusHalf)
    {
      LOOP_OVER_TILES(&ts->primNPlusHalf)
      {
        struct gridTile tile;
        setTile(iTile, jTile, kGlobal, &ts->primNPlusHalf, &tile);

        LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
        {
          struct gridZone zone;
          setZone(iInTile, jInTile, 0, &tile, &zone);

          REAL XCoords[NDIM];
          REAL primVars[DOF], sources[DOF];
          REAL connectionCoeffs[64];

          getXCoords(&zone, CENTER, XCoords);
          struct geometry geom; setGeometry(XCoords, &geom);

          memcpy(primVars, &INDEX_GRID(&ts->primNPlusHalf, &zone, 0),
                 sizeof(REAL[DOF])
                );

          struct fluidElement elem;
          setFluidElement(primVars, &geom, &elem);

          memcpy(connectionCoeffs, &INDEX_GRID(&ts->connection, &zone, 0),
                 sizeof(REAL[64])
                );

          computeSourceTerms(&elem, &geom, connectionCoeffs, sources);

          memcpy(&INDEX_GRID(&ts->sourcesNPlusHalf, &zone, 0), sources,
                 sizeof(REAL[DOF])
                );
        }
      }
    }

    finishFillingVecGhost(prim);

    setPointerToVecGhost(prim);

    LOOP_OVER_TILES(prim)
    {
      struct gridTile tile;
      setTile(iTile, jTile, kGlobal, &ts->primN, &tile);

      REAL primTile[TILE_SIZE(DOF)];
      REAL fluxX1Tile[TILE_SIZE(DOF)], fluxX2Tile[TILE_SIZE(DOF)];

      /* Load data from the global memory on RAM onto a tile small enough
        * to reside on the cache */
      loadPrimTile(prim, tile, primTile);

      /* Sync point */
      
//      applyTileBoundaryConditions(tile, primTile);
//
//      applyAdditionalProblemSpecificBCs(tile,
//                                        ts->problemSpecificData,
//                                        primTile
//                                       );
//      /* Sync point */
//  
//      /* Work on the tiles.*/
      computeFluxesOverTile(tile, primTile,
                            fluxX1Tile, fluxX2Tile,
                            ts->dtGrid
                           );

//      applyProblemSpecificFluxFilter(tile,
//                                     ts->problemSpecificData,
//                                     fluxX1Tile, fluxX2Tile
//                                    );
//
//      LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
//      {
//        struct gridZone zone;
//        setGridZone(iInTile, jInTile,
//                    &tile, &zone
//                   );
//
//        REAL XCoords[NDIM];
//        REAL primVars[DOF], sourceTerms[DOF], conservedVars[DOF];
//
//        getXCoords(&zone, CENTER, XCoords);
//        struct geometry geom; setGeometry(XCoords, &geom);
//
//        /* Now we need to compute conservedVarsOld using data from
//         * primOldGlobal. The computation of conservedVarsOld need not be done
//         * during the second half step. Put a switch here to avoid it. */
//        struct fluidElement elem;
//        setFluidElement(&INDEX_PETSC(primOldGlobal, &zone, 0), &geom, &elem);
//        computeFluxes(&elem, &geom, 0, conservedVars);
//
//        for (int var=0; var<DOF; var++)
//        {
//          INDEX_PETSC(conservedVarsOldGlobal, &zone, var) = 
//            conservedVars[var];
//
//          INDEX_PETSC(divFluxOldGlobal, &zone, var) = 
//            (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(&zone, var)]
//             - fluxX1Tile[INDEX_TILE(&zone, var)]
//            )/zone.dX1
//          #if (COMPUTE_DIM==2)
//          + 
//            (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(&zone, var)]
//             - fluxX2Tile[INDEX_TILE(&zone, var)]
//            )/zone.dX2
//          #endif
//            ;
//        }
//
//        if (ts->computeSourceTermsAtTimeN)
//        {
//          computeSourceTerms(&elem, &geom,
//                             &INDEX_PETSC(connectionGlobal, &zone, 0),
//                             sourceTerms);
//
//          for (int var=0; var<DOF; var++)
//          {
//            INDEX_PETSC(sourceTermsOldGlobal, &zone, var) = 
//              sourceTerms[var];
//          }
//        }
//        else if (ts->computeSourceTermsAtTimeNPlusHalf)
//        {
//          setFluidElement(&INDEX_PETSC(primHalfStepLocal, &zone, 0),
//                          &geom, &elem);
//
//          computeSourceTerms(&elem, &geom,
//                             &INDEX_PETSC(connectionGlobal, &zone, 0),
//                             sourceTerms);
//
//          for (int var=0; var<DOF; var++)
//          {
//            INDEX_PETSC(sourceTermsOldGlobal, &zone, var) = 
//              sourceTerms[var];
//          }
//        }
//
//      }
//
    } /* End of LOOP_OVER_TILES */

    restorePointerToVec(&ts->primN);
    restorePointerToVec(&ts->primNPlusHalf);

    /* All old sources and divFluxes have now been computed */
    ts->isZerothIterationOfSNES = 0;
  }
//
//  /* The following computation requires no communication*/
//
//  #if (TIME_STEPPING==IMPLICIT)
//    Vec primPetscVecLocal;
//    DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecLocal);
//
//    /* Exchange ghost zone data. */
//    DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
//                         primPetscVec,
//                         INSERT_VALUES,
//                         primPetscVecLocal);
//    DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
//                       primPetscVec,
//                       INSERT_VALUES,
//                       primPetscVecLocal);
//
//    ARRAY(primLocal);
//
//    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecLocal, &primLocal);
//  #endif
//
//  LOOP_OVER_TILES(X1Size, X2Size)
//  {
//    REAL primTile[TILE_SIZE];
//    REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];
//
//    #if (TIME_STEPPING==IMPLICIT)
//      LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
//      {
//        struct gridZone zone;
//        setGridZone(iTile, jTile,
//                    iInTile, jInTile,
//                    X1Start, X2Start, 
//                    X1Size, X2Size, 
//                    &zone);
//        for (int var=0; var<DOF; var++)
//        {
//          primTile[INDEX_TILE(&zone, var)] =
//          INDEX_PETSC(primLocal, &zone, var);
//        }
//      }
//      /* Sync point */
//    
//      /* Apply boundary conditions on each tile */
//      applyTileBoundaryConditions(iTile, jTile,
//                                  X1Start, X2Start,
//                                  X1Size, X2Size,
//                                  primTile);
//
//      applyAdditionalProblemSpecificBCs(iTile, jTile,
//                                        X1Start, X2Start,
//                                        X1Size, X2Size,
//                                        ts->problemSpecificData,
//                                        primTile);
//
//      computeFluxesOverTile(primTile, 
//                            iTile, jTile,
//                            X1Start, X2Start,
//                            X1Size, X2Size,
//                            fluxX1Tile, fluxX2Tile,
//                            dtGlobal);
//
//      applyProblemSpecificFluxFilter(iTile, jTile,
//                                     X1Start, X2Start,
//                                     X1Size, X2Size,
//                                     ts->problemSpecificData,
//                                     fluxX1Tile, fluxX2Tile);
//    #endif
//
//    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
//    {
//      struct gridZone zone;
//      setGridZone(iTile, jTile,
//                  iInTile, jInTile,
//                  X1Start, X2Start, 
//                  X1Size, X2Size, 
//                  &zone);
//
//      REAL XCoords[NDIM];
//
//      getXCoords(&zone, CENTER, XCoords);
//      struct geometry geom; setGeometry(XCoords, &geom);
//
//      struct fluidElement elem;
//      setFluidElement(&INDEX_PETSC(primGlobal, &zone, 0), &geom, &elem);
//
//      REAL conservedVars[DOF];
//      computeFluxes(&elem, &geom, 0, conservedVars);
//
//      #if (TIME_STEPPING==IMEX || TIME_STEPPING==IMPLICIT)
//        REAL sourceTerms[DOF];
//        computeSourceTerms(&elem, &geom,
//                           &INDEX_PETSC(connectionGlobal, &zone, 0),
//                           sourceTerms);
//      #endif
//
//      REAL g = sqrt(-geom.gDet);
//      REAL norm = g;
//
//      for (int var=0; var<DOF; var++)
//      {
//        #if (TIME_STEPPING==EXPLICIT)
//
//          INDEX_PETSC(residualGlobal, &zone, var) = 
//          ( (  conservedVars[var]
//             - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
//            )/ts->dt
//            + INDEX_PETSC(divFluxOldGlobal, &zone, var)
//            - INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
//          )/norm;
//	
//        #elif (TIME_STEPPING==IMEX)
//
//          INDEX_PETSC(residualGlobal, &zone, var) = 
//          ( (  conservedVars[var]
//             - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
//            )/ts->dt
//            + INDEX_PETSC(divFluxOldGlobal, &zone, var)
//            - 0.5*(  INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
//                   + sourceTerms[var]
//                  )
//	    )/norm;
//
//        #elif (TIME_STEPPING==IMPLICIT)
//		
//          INDEX_PETSC(residualGlobal, &zone, var) = 
//          ( (  conservedVars[var]
//             - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
//            )/ts->dt
//            + 0.5*(  INDEX_PETSC(divFluxOldGlobal, &zone, var)
//                   + 
//                    (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(&zone, var)]
//                     - fluxX1Tile[INDEX_TILE(&zone, var)]
//                    )/zone.dX1
//                  #if (COMPUTE_DIM==2)
//                   + 
//                    (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(&zone, var)]
//                     - fluxX2Tile[INDEX_TILE(&zone, var)]
//                    )/zone.dX2
//                  #endif
//                  )
//            - 0.5*(  INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
//                   + sourceTerms[var]
//                  )
//          )/norm;
//
//        #endif
//      }
//
//    }
//
//  } /* End of LOOP_OVER_TILES */
//
//
//
  return(0);
}
