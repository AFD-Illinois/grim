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

      INDEX_PETSC(primOldGlobal, &zone, RHO) = 25./(36.*M_PI);
      INDEX_PETSC(primOldGlobal, &zone, UU) = 5./(12.*M_PI*(ADIABATIC_INDEX - 1.));

      REAL v1 = -0.5*sin(2*M_PI*XCoords[2]);
      REAL v2 = 0.5*sin(2*M_PI*XCoords[1]);
      REAL v3 = 0.;
      REAL gamma = 1./sqrt(1 - v1*v1 - v2*v2 - v3*v3);

      INDEX_PETSC(primOldGlobal, &zone, U1) = gamma*v1;
      INDEX_PETSC(primOldGlobal, &zone, U2) = gamma*v2;
      INDEX_PETSC(primOldGlobal, &zone, U3) = gamma*v3;
      INDEX_PETSC(primOldGlobal, &zone, B1) = 
        -1./sqrt(4*M_PI)*sin(2*M_PI*XCoords[2]);
      INDEX_PETSC(primOldGlobal, &zone, B2) =
        1./sqrt(4*M_PI)*sin(4*M_PI*XCoords[1]);
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
  REAL kappaProblem=0.1;
  REAL tauProblem=.1; 

  void setConductionParameters(const struct geometry geom[ARRAY_ARGS 1],
                               struct fluidElement elem[ARRAY_ARGS 1])
  {
//    elem->kappa = kappaProblem;
//    elem->tau   = tauProblem;
  
    REAL xCoords[NDIM];
    XTox(geom->XCoords, xCoords);

    REAL P   = (ADIABATIC_INDEX-1.)*elem->primVars[UU];
    REAL T   = P/elem->primVars[RHO];
    REAL cs  = sqrt(  (ADIABATIC_INDEX-1.)*P
                    / (elem->primVars[RHO] + elem->primVars[UU])
                   );

    REAL phiCeil = elem->primVars[RHO] * pow(cs, 3.);

    REAL tauDynamical = 1.;
    REAL lambda       = 0.01;
    REAL y            = elem->primVars[PHI]/phiCeil;
    REAL fermiDirac   = 1./(exp((y-1.)/lambda) + 1.) + 1e-3;
    
    REAL tau    = tauDynamical*fermiDirac;
    elem->kappa = cs*cs*tau*elem->primVars[RHO];
    elem->tau   = tau; 
  }
#endif

#if (VISCOSITY)
  REAL etaProblem=0.1;
  REAL tauVisProblem=.1; 

  void setViscosityParameters(const struct geometry geom[ARRAY_ARGS 1],
                               struct fluidElement elem[ARRAY_ARGS 1])
  {
//    elem->eta     = etaProblem;
//    elem->tauVis  = tauVisProblem;

    REAL P   = (ADIABATIC_INDEX-1.)*elem->primVars[UU];
    REAL T   = P/elem->primVars[RHO];
    REAL cs  = sqrt(  (ADIABATIC_INDEX-1.)*P
                    / (elem->primVars[RHO] + elem->primVars[UU])
                   );

    REAL bSqr    = getbSqr(elem, geom);
    REAL beta    = P/(bSqr/2.);
    REAL psiCeil = 2.*P/beta;

    REAL tauDynamical = 1.;
    REAL lambda       = 0.01;
    REAL y            = fabs(elem->primVars[PSI])/fabs(psiCeil);
    REAL fermiDirac   = 1./(exp((y-1.)/lambda) + 1.) + 1e-3;
    
    REAL tau     = tauDynamical*fermiDirac;
    elem->eta    = cs*cs*tau*elem->primVars[RHO];
    elem->tauVis = tau;
  }
#endif

void writeProblemSpecificData(PetscViewer parametersViewer,
  const struct problemData problemSpecificData[ARRAY_ARGS 1]) {

}
