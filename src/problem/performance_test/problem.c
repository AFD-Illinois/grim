#include "../problem.h"
#include "torus.h"

#if (CONDUCTION || VISCOSITY)
  void setDiffusionCoefficients(const struct geometry geom[ARRAY_ARGS 1],
                                struct fluidElement elem[ARRAY_ARGS 1]
                               )
  {
    #if (VISCOSITY)
      elem->nu  = 0.1;
    #endif
    #if (CONDUCTION)
      elem->chi = 0.1; 
    #endif
  }
#endif

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  PetscRandom randomNumberGenerator;
  PetscRandomCreate(PETSC_COMM_WORLD, &randomNumberGenerator);
  VecSetRandom(ts->primN.vec, randomNumberGenerator);
  PetscRandomDestroy(&randomNumberGenerator);

  PetscPrintf(PETSC_COMM_WORLD, "In here\n");
}

//void applyFloor(const int iTile, const int jTile,
//                const int X1Start, const int X2Start,
//                const int X1Size, const int X2Size,
//                REAL primTile[ARRAY_ARGS TILE_SIZE])
//{
//}
//
//void applyAdditionalProblemSpecificBCs
//(
//  const int iTile, const int jTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  const struct problemData problemSpecificData[ARRAY_ARGS 1],
//  REAL primTile[ARRAY_ARGS TILE_SIZE]
//)
//{
//}
//
//void applyProblemSpecificFluxFilter
//(
//  const int iTile, const int jTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  const struct problemData problemSpecificData[ARRAY_ARGS 1],
//  REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
//  REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]
//)
//{
//}

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
}

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
}
