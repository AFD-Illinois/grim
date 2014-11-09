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

      INDEX_PETSC(primOldGlobal, &zone, B1) = 0.;
      INDEX_PETSC(primOldGlobal, &zone, B2) = 0.;
      INDEX_PETSC(primOldGlobal, &zone, B3) = 0.;

      if (r < R)
      {
        INDEX_PETSC(primOldGlobal, &zone, B1) = -A0*XCoords[2]/r;

        INDEX_PETSC(primOldGlobal, &zone, B2) = A0*XCoords[1]/r;
      }

//      INDEX_PETSC(primOldGlobal, &zone, RHO) = 1 + exp(-r*r/0.01);
//      INDEX_PETSC(primOldGlobal, &zone, UU) = 1./(ADIABATIC_INDEX - 1.);
//      INDEX_PETSC(primOldGlobal, &zone, U1) = 4.95;
//      INDEX_PETSC(primOldGlobal, &zone, U2) = 4.95;
//      INDEX_PETSC(primOldGlobal, &zone, U3) = 0.;
//      INDEX_PETSC(primOldGlobal, &zone, B1) = 0.;
//      INDEX_PETSC(primOldGlobal, &zone, B2) = 0.;
//      INDEX_PETSC(primOldGlobal, &zone, B3) = 0.;

//      REAL M_PI = 3.14159265359;
//      INDEX_PETSC(primOldGlobal, &zone, RHO) = 25./(36.*M_PI);
//      INDEX_PETSC(primOldGlobal, &zone, UU) = 5./(12.*M_PI*(ADIABATIC_INDEX - 1.));
//
//      REAL v1 = -0.5*sin(2*M_PI*XCoords[2]);
//      REAL v2 = 0.5*sin(2*M_PI*XCoords[1]);
//      REAL v3 = 0.;
//      REAL gamma = 1./sqrt(1 - v1*v1 - v2*v2 - v3*v3);
//
//      INDEX_PETSC(primOldGlobal, &zone, U1) = gamma*v1;
//      INDEX_PETSC(primOldGlobal, &zone, U2) = gamma*v2;
//      INDEX_PETSC(primOldGlobal, &zone, U3) = gamma*v3;
//      INDEX_PETSC(primOldGlobal, &zone, B1) = 
//        -1./sqrt(4*M_PI)*sin(2*M_PI*XCoords[2]);
//      INDEX_PETSC(primOldGlobal, &zone, B2) =
//        1./sqrt(4*M_PI)*sin(4*M_PI*XCoords[1]);
//      INDEX_PETSC(primOldGlobal, &zone, B3) = 0.;
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
