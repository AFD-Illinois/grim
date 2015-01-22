#include "../problem.h"

#if (CONDUCTION)
  REAL kappaProblem;
  REAL tauProblem; 

  void setConductionParameters(struct fluidElement elem[ARRAY_ARGS 1])
  {
    elem->kappa = kappaProblem;
    elem->tau   = tauProblem;
  }
#endif

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecOld, &primOldGlobal);

  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL XCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);


      #if (MODE==HYDRO_SOUND_MODE_1D) /* Eigenvalue = 3.09362659024*I */
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;

        deltaPrimVars[RHO] = 0.345991032308;
        deltaPrimVars[UU]  = 0.922642752822;
        deltaPrimVars[U1]  = -0.170354208129;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }

      #elif (MODE==ENTROPY_WAVE_1D) /* Eigenvalue = 0 */
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;

        deltaPrimVars[RHO] = 1.;
        deltaPrimVars[UU]  = 0.;
        deltaPrimVars[U1]  = 0.;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }

      #elif (MODE==CONDUCTION_STABLE_1D) 
        /* Eigenvalue = -0.498597173331 - 0.857832357798*I
         * kappa = 0.1
         * tau = 1.01818181818182 */
        kappaProblem = .1;
        tauProblem   = 1.01818181818182;

        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.000001;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
        primVars0[PHI] = 0.;

        deltaPrimVars[RHO] = 0.911376933183;
        deltaPrimVars[UU]  = 0.030751595371 - 0.0635975709194*I;
        deltaPrimVars[U1]  = 0.124428706971 - 0.0723215917578*I;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PHI] = -0.332658158111 + 0.181734443922*I;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #endif

    }
  }
  
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld, &primOldGlobal);
}


void applyFloor(const int iTile, const int jTile,
                const int X1Start, const int X2Start,
                const int X1Size, const int X2Size,
                REAL primTile[ARRAY_ARGS TILE_SIZE])
{
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    if (primTile[INDEX_TILE(&zone, RHO)] < RHO_FLOOR)
    {
        primTile[INDEX_TILE(&zone, RHO)] = RHO_FLOOR;
    }

    if (primTile[INDEX_TILE(&zone, UU)] < UU_FLOOR)
    {
        primTile[INDEX_TILE(&zone, UU)] = UU_FLOOR;
    }
  }
}

void applyFloorInsideNonLinearSolver(const int iTile, const int jTile,
                                     const int X1Start, const int X2Start,
                                     const int X1Size, const int X2Size,
                                     REAL primTile[ARRAY_ARGS TILE_SIZE])
{

}

void applyAdditionalProblemSpecificBCs(const int iTile, const int jTile,
                                       const int X1Start, const int X2Start,
                                       const int X1Size, const int X2Size,
                                       REAL tile[ARRAY_ARGS TILE_SIZE])
{

}

void applyProblemSpecificFluxFilter(const int iTile, const int jTile,
                                    const int X1Start, const int X2Start,
                                    const int X1Size, const int X2Size,
                                    REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
                                    REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE])
{

}

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{

}

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{

}

