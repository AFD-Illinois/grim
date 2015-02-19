#include "../problem.h"

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

      REAL X1Center;
      X1Center = (X1_A + X1_B)/2.;

      REAL rhoLeft,       rhoRight;
      REAL pressureLeft,  pressureRight;
      REAL u1Left,        u1Right;
      REAL u2Left,        u2Right;
      REAL u3Left,        u3Right;
      REAL B1Left,        B1Right;
      REAL B2Left,        B2Right;
      REAL B3Left,        B3Right;

      #if (SHOCK_TEST == FAST_SHOCK)

        rhoLeft       = 1.;         rhoRight      = 25.48;
        pressureLeft  = 1.;         pressureRight = 367.5;
        u1Left        = 25.;        u1Right       = 1.091;
        u2Left        = 0.;         u2Right       = 0.3923;
        u3Left        = 0.;         u3Right       = 0.;
        B1Left        = 20.;        B1Right       = 20.;
        B2Left        = 25.02;      B2Right       = 49.;
        B3Left        = 0.;         B3Right       = 0.;

      #elif (SHOCK_TEST == SLOW_SHOCK)

        rhoLeft       = 1.;         rhoRight      = 3.323;
        pressureLeft  = 10.;        pressureRight = 55.36;
        u1Left        = 1.53;       u1Right       = 0.9571;
        u2Left        = 0.;         u2Right       = -0.6822;
        u3Left        = 0.;         u3Right       = 0.;
        B1Left        = 10.;        B1Right       = 10.;
        B2Left        = 18.28;      B2Right       = 14.49;
        B3Left        = 0.;         B3Right       = 0.;

      #elif (SHOCK_TEST == SWITCH_ON_SLOW)

        rhoLeft       = 1.78e-3;    rhoRight      = 0.01;
        pressureLeft  = .1;         pressureRight = 1.;
        u1Left        = -.765;      u1Right       = 0.;
        u2Left        = -1.386;     u2Right       = 0.;
        u3Left        = 0.;         u3Right       = 0.;
        B1Left        = 1.;         B1Right       = 1.;
        B2Left        = 1.022;      B2Right       = 0.;
        B3Left        = 0.;         B3Right       = 0.;

      #elif (SHOCK_TEST == SWITCH_OFF_FAST)

        rhoLeft       = 0.1;        rhoRight      = 0.562;
        pressureLeft  = 1.;         pressureRight = 10.;
        u1Left        = -2.;        u1Right       = -0.212;
        u2Left        = 0.;         u2Right       = -0.590;
        u3Left        = 0.;         u3Right       = 0.;
        B1Left        = 2.;         B1Right       = 2.;
        B2Left        = 0.;         B2Right       = 4.71;
        B3Left        = 0.;         B3Right       = 0.;

      #elif (SHOCK_TEST == ALFVEN_WAVE)

        rhoLeft       = 1.;         rhoRight      = 1.;
        pressureLeft  = 1.;         pressureRight = 1.;
        u1Left        = 0.;         u1Right       = 3.7;
        u2Left        = 0.;         u2Right       = 5.76;
        u3Left        = 0.;         u3Right       = 0.;
        B1Left        = 3.;         B1Right       = 3.;
        B2Left        = 3.;         B2Right       = -6.857;
        B3Left        = 0.;         B3Right       = 0.;

      #elif (SHOCK_TEST == SHOCK_TUBE_1)

        rhoLeft       = 1.;         rhoRight      = .1;
        pressureLeft  = 1000.;      pressureRight = 1.;
        u1Left        = 0.;         u1Right       = 0.;
        u2Left        = 0.;         u2Right       = 0.;
        u3Left        = 0.;         u3Right       = 0.;
        B1Left        = 1.;         B1Right       = 1.;
        B2Left        = 0.;         B2Right       = 0.;
        B3Left        = 0.;         B3Right       = 0.;

      #elif (SHOCK_TEST == SHOCK_TUBE_2)

        rhoLeft       = 1.;         rhoRight      = .1;
        pressureLeft  = 30.;        pressureRight = 1.;
        u1Left        = 0.;         u1Right       = 0.;
        u2Left        = 0.;         u2Right       = 0.;
        u3Left        = 0.;         u3Right       = 0.;
        B1Left        = 0.;         B1Right       = 0.;
        B2Left        = 20.;        B2Right       = 0.;
        B3Left        = 0.;         B3Right       = 0.;

      #elif (SHOCK_TEST == COLLISION)
      
        rhoLeft       = 1.;         rhoRight      = 1.;
        pressureLeft  = 1.;         pressureRight = 1.;
        u1Left        = 5.;         u1Right       = -5.;
        u2Left        = 0.;         u2Right       = 0.;
        u3Left        = 0.;         u3Right       = 0.;
        B1Left        = 10.;        B1Right       = 10.;
        B2Left        = 10.;        B2Right       = -10.;
        B3Left        = 0.;         B3Right       = 0.;

      #elif (SHOCK_TEST == STATIONARY_SHOCK)
      
        rhoLeft       = 1.;                      rhoRight      = 3.08312999;
        pressureLeft  = 1./(ADIABATIC_INDEX-1.); pressureRight = 4.94577705/(ADIABATIC_INDEX-1.);
        u1Left        = 1.;         u1Right       = 0.32434571;
        u2Left        = 0.;         u2Right       = 0.;
        u3Left        = 0.;         u3Right       = 0.;
        B1Left        = 0.00001;    B1Right       = 0.00001;
        B2Left        = 0.;         B2Right       = 0.;
        B3Left        = 0.;         B3Right       = 0.;


      #endif

      if (XCoords[1] < X1Center)
      {
        INDEX_PETSC(primOldGlobal, &zone, RHO)  = rhoLeft;
        INDEX_PETSC(primOldGlobal, &zone, UU)   = pressureLeft/(ADIABATIC_INDEX-1.);
        INDEX_PETSC(primOldGlobal, &zone, U1)   = u1Left;
        INDEX_PETSC(primOldGlobal, &zone, U2)   = u2Left;
        INDEX_PETSC(primOldGlobal, &zone, U3)   = u3Left;
        INDEX_PETSC(primOldGlobal, &zone, B1)   = B1Left;
        INDEX_PETSC(primOldGlobal, &zone, B2)   = B2Left;
        INDEX_PETSC(primOldGlobal, &zone, B3)   = B3Left;
      } else
      {
        INDEX_PETSC(primOldGlobal, &zone, RHO)  = rhoRight;
        INDEX_PETSC(primOldGlobal, &zone, UU)   = pressureRight/(ADIABATIC_INDEX-1.);
        INDEX_PETSC(primOldGlobal, &zone, U1)   = u1Right;
        INDEX_PETSC(primOldGlobal, &zone, U2)   = u2Right;
        INDEX_PETSC(primOldGlobal, &zone, U3)   = u3Right;
        INDEX_PETSC(primOldGlobal, &zone, B1)   = B1Right;
        INDEX_PETSC(primOldGlobal, &zone, B2)   = B2Right;
        INDEX_PETSC(primOldGlobal, &zone, B3)   = B3Right;
      }

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
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    if (primTile[INDEX_TILE(&zone, RHO)] < RHO_FLOOR_MIN)
    {
      primTile[INDEX_TILE(&zone, RHO)] = RHO_FLOOR_MIN;
    }

    if (primTile[INDEX_TILE(&zone, UU)] < UU_FLOOR_MIN)
    {
      primTile[INDEX_TILE(&zone, UU)] = UU_FLOOR_MIN;
    }
  }
}

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primHalfStepGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecHalfStep, &primHalfStepGlobal);

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] =
        INDEX_PETSC(primHalfStepGlobal, &zone, var);
      }
    }

    applyFloor(iTile, jTile,
               ts->X1Start, ts->X2Start,
               ts->X1Size, ts->X2Size,
               primTile);

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primHalfStepGlobal, &zone, var) =
        primTile[INDEX_TILE(&zone, var)];
      }
    }

  }

  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecHalfStep, &primHalfStepGlobal);
}

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecOld, &primOldGlobal);

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] =
        INDEX_PETSC(primOldGlobal, &zone, var);
      }
    }

    applyFloor(iTile, jTile,
               ts->X1Start, ts->X2Start,
               ts->X1Size, ts->X2Size,
               primTile);

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primOldGlobal, &zone, var) =
        primTile[INDEX_TILE(&zone, var)];
      }
    }

  }

  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld, &primOldGlobal);
}

void applyAdditionalProblemSpecificBCs
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL primTile[ARRAY_ARGS TILE_SIZE]
)
{

}

void applyProblemSpecificFluxFilter
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
  REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]
)
{

}

#if (CONDUCTION)
void setConductionParameters(const struct geometry geom[ARRAY_ARGS 1],
                             struct fluidElement elem[ARRAY_ARGS 1])
{
  SETERRQ(PETSC_COMM_WORLD, 1,
          "Conduction parameters not set in shock_tests/problem.c\n");
}
#endif
