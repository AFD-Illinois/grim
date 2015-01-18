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
        /* Eigenvalue = -4.99386690579 - 13.32147544*I
         * kappa = 2.
         * tau = 0.463636363636364 */
        kappaProblem = 2.;
        tauProblem   = 0.463636363636364;

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
        primVars0[PHI] = 0.;

        deltaPrimVars[RHO] = -0.106491811933 - 0.0383867866432*I;
        deltaPrimVars[UU]  = 0.109942056865 + 0.0453054519637*I;
        deltaPrimVars[U1]  = -0.256291432121 + 0.0032526972965*I;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.00001;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PHI] = 0.952549332339;

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

