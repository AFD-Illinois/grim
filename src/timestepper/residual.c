#include "timestepper.h"

PetscErrorCode computeResidual(SNES snes, 
                               Vec primVec,
                               Vec residualVec,
                               void *ptr
                              )
{
  struct timeStepper *ts = (struct timeStepper*)ptr;
  struct gridData *prim  = &ts->primN;

  setPointerToExternalVec(primVec, &ts->primNPlusOne);
  setPointerToVec(&ts->primN);
  setPointerToVec(&ts->primNPlusHalf);
  setPointerToVec(&ts->conservedVarsN);
  setPointerToVec(&ts->sources);
  setPointerToVec(&ts->divFluxes);
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

    if (ts->computeConsVarsAndSourcesAtN)
    {  
      LOOP_OVER_TILES(&ts->primN)
      {
        struct gridTile tile;
        SET_TILE_INDICES(iTile, jTile, kGlobal, &tile);

        LOOP_INSIDE_TILE(0, TILE_SIZE_X2, 0, TILE_SIZE_X3)
        {
          struct gridStrip strip;
          SET_STRIP_INDICES(jInTile, kInTile, &strip);

          struct fluidElement elem;
          setFluidElement(prim, &strip, geom, &elem);

          computeConservedVars(&elem, geom, conservedVarsN);
          computeSources(&elem, geom, sources);
        }
      }
    }

    if (ts->computeSourcesAtNPlusHalf)
    {
      LOOP_OVER_TILES(&ts->primNPlusHalf)
      {
        struct gridTile tile;
        SET_TILE_INDICES(iTile, jTile, kGlobal, &tile);

        LOOP_INSIDE_TILE(0, TILE_SIZE_X2, 0, TILE_SIZE_X3)
        {
          struct gridStrip strip;
          SET_STRIP_INDICES(jInTile, kInTile, &strip);

          struct fluidElement elem;
          setFluidElement(primNPlusHalf, &strip, geom, &elem);

          computeSources(&elem, geom, sources);
        }
      }
    }

    finishFillingVecGhost(prim);

    setPointerToVecGhost(prim);

    LOOP_OVER_TILES(prim)
    {
      struct gridTile tile;
      setTile(iTile, jTile, kGlobal, prim, &tile);

      REAL primTile[TILE_SIZE(DOF)];
      REAL fluxesTile[COMPUTE_DIM][TILE_SIZE(DOF)];

      /* Load data from the global memory on RAM onto a tile small enough
        * to reside on the cache */
      loadPrimTile(prim, &tile, primTile);

      /* Sync point */
      
      applyTileBoundaryConditions(&tile, primTile);

//      applyAdditionalProblemSpecificBCs(tile,
//                                        ts->problemSpecificData,
//                                        primTile
//                                       );
      /* Sync point */
  
      /* Work on the tiles.*/
      computeFluxesOverTile(&tile, primTile, fluxesTile,
                            &ts->dtGrid
                           );

//      applyProblemSpecificFluxFilter(tile,
//                                     ts->problemSpecificData,
//                                     fluxX1Tile, fluxX2Tile
//                                    );

      LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2, 0, TILE_SIZE_X3)
      {
        struct gridZone zone;
        setZone(iInTile, jInTile, kInTile, &tile, &zone);

        REAL divFluxes[DOF];

        for (int var=0; var<DOF; var++)
        {
          divFluxes[var] =
            (  fluxesTile[X1][INDEX_TILE_OFFSET(1, 0, 0, &zone, var)]
             - fluxesTile[X1][INDEX_TILE(&zone, var)]
            )/zone.dX[1]
          #if (COMPUTE_DIM==2)
          + 
            (  fluxesTile[X2][INDEX_TILE_OFFSET(0, 1, 0, &zone, var)]
             - fluxesTile[X2][INDEX_TILE(&zone, var)]
            )/zone.dX[2]
          #endif
          #if (COMPUTE_DIM==3)
          + 
            (  fluxesTile[X3][INDEX_TILE_OFFSET(0, 0, 1, &zone, var)]
             - fluxesTile[X3][INDEX_TILE(&zone, var)]
            )/zone.dX[3]
          #endif
            ;
        }

        memcpy(&INDEX_GRID(&ts->divFluxes, &zone, 0), divFluxes,
               sizeof(REAL[DOF])
              );
      }

    } /* End of LOOP_OVER_TILES */

    restorePointerToVecGhost(prim);

    /* All old sources and divFluxes have now been computed */
    ts->isZerothIterationOfSNES = 0;
  }

  /* The following computation requires no communication except if using
   * TIME_STEPPING==IMPLICIT */

  #if (TIME_STEPPING==IMPLICIT)
    startFillingVecGhostWithExternalVec(primVec, &ts->primNPlusOne);
  #endif

  LOOP_OVER_TILES(&ts->residual)
  {
    struct gridTile tile;
    setTile(iTile, jTile, kGlobal, &ts->residual, &tile);

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2, 0, TILE_SIZE_X3)
    {
      struct gridZone zone;
      setZone(iInTile, jInTile, kInTile, &tile, &zone);

      REAL XCoords[NDIM];
      REAL primVarsNPlusOne[DOF];
      REAL conservedVarsN[DOF], conservedVarsNPlusOne[DOF];
      REAL sources[DOF], divFluxes[DOF];
      REAL residual[DOF];

      memcpy(primVarsNPlusOne, &INDEX_GRID(&ts->primNPlusOne, &zone, 0),
             sizeof(REAL[DOF])
            );
      memcpy(conservedVarsN, &INDEX_GRID(&ts->conservedVarsN, &zone, 0),
             sizeof(REAL[DOF])
            );
      memcpy(divFluxes, &INDEX_GRID(&ts->divFluxes, &zone, 0),
             sizeof(REAL[DOF])
            );
      memcpy(sources, &INDEX_GRID(&ts->sources, &zone, 0),
             sizeof(REAL[DOF])
            );

      getXCoords(&zone, CENTER, XCoords);
      struct geometry geom; setGeometry(XCoords, &geom);

      struct fluidElement elem;
      setFluidElement(primVarsNPlusOne, &geom, &elem);
      computeFluxes(&elem, &geom, 0, conservedVarsNPlusOne);

      #if (TIME_STEPPING == IMPLICIT || TIME_STEPPING == IMEX)
        REAL sourcesNPlusOne[DOF];
        REAL connectionCoeffs[64];
        memcpy(connectionCoeffs, &INDEX_GRID(&ts->connection, &zone, 0),
               sizeof(REAL[64])
              );
        computeSourceTerms(&elem, &geom, connectionCoeffs, sourcesNPlusOne);
      #endif

      REAL g = sqrt(-geom.gDet);
      REAL norm = g;

      for (int var=0; var<DOF; var++)
      {
        residual[var] =  
          (  conservedVarsNPlusOne[var] 
           - conservedVarsN[var]
          )/ts->dt
               
        #if (TIME_STEPPING == EXPLICIT)

          + divFluxes[var] - sources[var]

        #elif (TIME_STEPPING == IMEX)
                        
          + divFluxes[var] - 0.5*(sources[var] + sourcesNPlusOne[var])

        #elif (TIME_STEPPING == IMPLICIT)
        
          +  0.5*divFluxes[var] - 0.5*(sources[var] + sourcesNPlusOne[var])
        /* We will add in the remaining 0.5*divFluxes[var] for IMPLICIT below */

        #endif
        ;
      }

      memcpy(&INDEX_GRID(&ts->residual, &zone, 0), residual, 
             sizeof(REAL[DOF])
            );
    }
  }

  #if (TIME_STEPPING==IMPLICIT)
    finishFillingVecGhostWithExternalVec(primVec, &ts->primNPlusOne);
    setPointerToVecGhost(&ts->primNPlusOne);

    LOOP_OVER_TILES(&ts->primNPlusOne)
    {
      struct gridTile tile;
      setTile(iTile, jTile, kGlobal, &ts->primN, &tile);

      REAL primTile[TILE_SIZE(DOF)];
      REAL fluxesTile[COMPUTE_DIM][TILE_SIZE(DOF)];

      /* Load data from the global memory on RAM onto a tile small enough
        * to reside on the cache */
      loadPrimTile(&ts->primNPlusOne, &tile, primTile);

      /* Sync point */
      
      applyTileBoundaryConditions(&tile, primTile);

//      applyAdditionalProblemSpecificBCs(tile,
//                                        ts->problemSpecificData,
//                                        primTile
//                                       );
      /* Sync point */
  
      /* Work on the tiles.*/
      computeFluxesOverTile(&tile, primTile, fluxesTile,
                            &ts->dtGrid
                           );

//      applyProblemSpecificFluxFilter(tile,
//                                     ts->problemSpecificData,
//                                     fluxX1Tile, fluxX2Tile
//                                    );

      LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2, 0, TILE_SIZE_X3)
      {
        struct gridZone zone;
        setZone(iInTile, jInTile, kInTile, &tile, &zone);

        REAL residual[DOF];
        memcpy(residual, &INDEX_GRID(&ts->residual, &zone, 0),
               sizeof(REAL[DOF])
              );

        for (int var=0; var<DOF; var++)
        {
          residual[var] +=
            0.5 * (
                    (  fluxesTile[X1][INDEX_TILE_OFFSET(1, 0, 0, &zone, var)]
                     - fluxesTile[X1][INDEX_TILE(&zone, var)]
                    )/zone.dX[1]
                  #if (COMPUTE_DIM==2)
                  + 
                    (  fluxesTile[X2][INDEX_TILE_OFFSET(0, 1, 0, &zone, var)]
                     - fluxesTile[X2][INDEX_TILE(&zone, var)]
                    )/zone.dX[2]
                  #endif
                  #if (COMPUTE_DIM==3)
                  + 
                    (  fluxesTile[X3][INDEX_TILE_OFFSET(0, 0, 1, &zone, var)]
                     - fluxesTile[X3][INDEX_TILE(&zone, var)]
                    )/zone.dX[3]
                  #endif
                );
        }

        memcpy(&INDEX_GRID(&ts->residual, &zone, 0), residual,
               sizeof(REAL[DOF])
              );
      }
    } /* End of LOOP_OVER_TILES */

    restorePointerToVecGhost(&ts->primNPlusOne);
  #endif

  restorePointerToExternalVec(primVec, &ts->primNPlusOne);
  restorePointerToVec(&ts->primN);
  restorePointerToVec(&ts->primNPlusHalf);
  restorePointerToVec(&ts->conservedVarsN);
  restorePointerToVec(&ts->sources);
  restorePointerToVec(&ts->divFluxes);
  restorePointerToVec(&ts->connection);
  restorePointerToVec(&ts->dtGrid);

  return(0);
}

void computeSourcesAndConservedVarsOverGrid
  (const int computeConservedVars,
   const struct gridData prim[ARRAY_ARGS 1],
   const struct gridData connection[ARRAY_ARGS 1],
   struct gridData sourcesGrid[ARRAY_ARGS 1],
   struct gridData conservedVarsGrid[ARRAY_ARGS 1]
  )
{
  if (computeConservedVars)
  {
    LOOP_OVER_TILES(prim)
    {
      struct gridTile tile;
      setTile(iTile, jTile, kGlobal, prim, &tile);

      LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2, 0, TILE_SIZE_X3)
      {
        struct gridZone zone;
        setZone(iInTile, jInTile, kInTile, &tile, &zone);

        REAL XCoords[NDIM];
        REAL primVars[DOF], sources[DOF], conservedVars[DOF];
        REAL connectionCoeffs[64];

        getXCoords(&zone, CENTER, XCoords);
        struct geometry geom; setGeometry(XCoords, &geom);

        /* Fast copy: Destination, Source, Size */
        memcpy(primVars, &INDEX_GRID(prim, &zone, 0), sizeof(REAL[DOF]) );

        struct fluidElement elem;
        setFluidElement(primVars, &geom, &elem);
        computeFluxes(&elem, &geom, 0, conservedVars);

        memcpy(connectionCoeffs, &INDEX_GRID(connection, &zone, 0),
               sizeof(REAL[64])
              );

        computeSourceTerms(&elem, &geom, connectionCoeffs, sources);

        memcpy(&INDEX_GRID(conservedVarsGrid, &zone, 0), conservedVars,
               sizeof(REAL[DOF])
              );

        memcpy(&INDEX_GRID(sourcesGrid, &zone, 0), sources,
               sizeof(REAL[DOF])
              );
      }
    }
  }
  else
  {
    LOOP_OVER_TILES(prim)
    {
      struct gridTile tile;
      setTile(iTile, jTile, kGlobal, prim, &tile);

      LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2, 0, TILE_SIZE_X3)
      {
        struct gridZone zone;
        setZone(iInTile, jInTile, kInTile, &tile, &zone);

        REAL XCoords[NDIM];
        REAL primVars[DOF], sources[DOF];
        REAL connectionCoeffs[64];

        getXCoords(&zone, CENTER, XCoords);
        struct geometry geom; setGeometry(XCoords, &geom);

        memcpy(primVars, &INDEX_GRID(prim, &zone, 0), sizeof(REAL[DOF]) );

        struct fluidElement elem;
        setFluidElement(primVars, &geom, &elem);

        memcpy(connectionCoeffs, &INDEX_GRID(connection, &zone, 0), 
               sizeof(REAL[64])
              );

        computeSourceTerms(&elem, &geom, connectionCoeffs, sources);

        memcpy(&INDEX_GRID(sourcesGrid, &zone, 0), sources,
               sizeof(REAL[DOF])
              );
      }
    }
  }
}
