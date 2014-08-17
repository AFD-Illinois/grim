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

      REAL X1Center, X2Center;
      X1Center = (X1_A + X1_B)/2.; X2Center = (X2_A + X2_B)/2.;

      REAL r = sqrt( pow((XCoords[1] - X1Center), 2.) +
                     pow((XCoords[2] - X2Center), 2.)
                   );

      INDEX_PETSC(primOldGlobal, &zone, RHO) = 1.;
      INDEX_PETSC(primOldGlobal, &zone, UU) = 1./(ADIABATIC_INDEX - 1.);

      REAL gamma = 1./sqrt(1. - V1*V1 - V2*V2 - V3*V3);

      INDEX_PETSC(primOldGlobal, &zone, U1) = gamma*V1;
      INDEX_PETSC(primOldGlobal, &zone, U2) = gamma*V2;
      INDEX_PETSC(primOldGlobal, &zone, U3) = gamma*V3;

      INDEX_PETSC(primOldGlobal, &zone, B1) =
        -2*A0*XCoords[2]/(R*R)*exp(-r*r/(R*R));

      INDEX_PETSC(primOldGlobal, &zone, B2) =
        -2*A0*XCoords[1]/(R*R)*exp(-r*r/(R*R));

      INDEX_PETSC(primOldGlobal, &zone, B3) = 0.;

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

void problemDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{

}
