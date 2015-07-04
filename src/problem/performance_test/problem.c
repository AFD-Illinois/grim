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
  VecSetRandom(ts->primN.vec, NULL);

//  startFillingVecGhost(&ts->primN);
//
//  finishFillingVecGhost(&ts->primN);
  setPointerToVec(&ts->primN);
  setPointerToVec(&ts->connection);
  setPointerToVec(&ts->conservedVarsN);
  setPointerToVec(&ts->sources);

  for (int n=0; n<1000; n++)
  {
    PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);
    computeSourcesAndConservedVarsOverGrid
      (1, &ts->primN, &ts->connection, &ts->sources, &ts->conservedVarsN);
  }

  restorePointerToVec(&ts->primN);
  restorePointerToVec(&ts->connection);
  restorePointerToVec(&ts->conservedVarsN);
  restorePointerToVec(&ts->sources);
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
