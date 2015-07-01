#ifndef GRIM_TIMESTEPPER_H_
#define GRIM_TIMESTEPPER_H_

#include <petsc.h>
#include <petscviewerhdf5.h>
#include <unistd.h>
#include "../inputs.h"
#include "../grid/grid.h"
#include "../geometry/geometry.h"
#include "../boundary/boundary.h"
#include "../reconstruct/reconstruct.h"
#include "../riemannsolver/riemannsolver.h"
#include "../physics/physics.h"
#include "../problem/problem.h"

#if (USE_OPENMP)
  #include <omp.h>
#endif

struct problemData;

struct timeStepper
{
  REAL t, dt, tDump;
  int timeStepCounter;
  int dumpCounter;

  SNES snes;

  struct gridData primNPlusOne;
  struct gridData primNPlusHalf;
  struct gridData primN;

  struct gridData conservedVarsN;
  struct gridData sourcesN;
  struct gridData sourcesNPlusHalf;
  
  struct gridData residual;

  struct gridData dtGrid;
  struct gridData connection;

  struct problemData *problemSpecificData;

  int isZerothIterationOfSNES;

  int computeDivOfFluxAtN;
  int computeDivOfFluxAtNPlusHalf;
  int computeSourcesAtN;
  int computeSourcesAtNPlusHalf;
  int computeConsVarsAndSourcesAtN;
};

/* User functions */

void timeStepperInit(struct timeStepper ts[ARRAY_ARGS 1]);

void timeStep(struct timeStepper ts[ARRAY_ARGS 1]);

void timeStepperDestroy(struct timeStepper ts[ARRAY_ARGS 1]);

/* Internal functions */

PetscErrorCode computeResidual(SNES snes, 
                               Vec primVec,
                               Vec residalVec,
                               void *ptr
                              );

//void computeFluxesOverTile(const REAL primTile[ARRAY_ARGS TILE_SIZE],
//                           const struct gridTile tile[ARRAY_ARGS 1],
//                           REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
//                           REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE],
//                           ARRAY(dtGlobal)
//                          );
//
//void setChristoffelSymbols(struct timeStepper ts[ARRAY_ARGS 1]);
//
//void diagnostics(struct timeStepper ts[ARRAY_ARGS 1]);
//
////#if (CONDUCTION || VISCOSITY)
////  void initDissipationDataStructures(struct timeStepper ts[ARRAY_ARGS 1]);
////  void destroyDissipationDataStructures(struct timeStepper ts[ARRAY_ARGS 1]);
////#endif
//
#endif /* GRIM_TIMESTEPPER_H_ */
