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
#include "../physics/physics.h"
#include "../problem/problem.h"
#include "macros.h"

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
  struct gridData divFluxes;  /* Stores either divFluxesN or divFluxesNPlusHalf */
  struct gridData sources;    /* Stores either sourcesN or sourcesNPlusHalf */
  
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

void computeFluxesOverTile
  (const struct gridTile tile[ARRAY_ARGS 1],
   const REAL primTile[ARRAY_ARGS TILE_SIZE(DOF)],
   REAL fluxesTile[ARRAY_ARGS COMPUTE_DIM][TILE_SIZE(DOF)],
   struct gridData dtGrid[ARRAY_ARGS 1]
  );

void fluxCT(const struct gridTile tile[ARRAY_ARGS 1],
            REAL fluxesTile[ARRAY_ARGS COMPUTE_DIM][TILE_SIZE(DOF)]
           );

void computeSourcesAndConservedVarsOverGrid
  (const int computeConservedVars,
   const struct gridData prim[ARRAY_ARGS 1],
   const struct gridData connection[ARRAY_ARGS 1],
   struct gridData sourcesGrid[ARRAY_ARGS 1],
   struct gridData conservedVarsGrid[ARRAY_ARGS 1]
  );

void setChristoffelSymbols(struct timeStepper ts[ARRAY_ARGS 1]);

//void diagnostics(struct timeStepper ts[ARRAY_ARGS 1]);
//
////#if (CONDUCTION || VISCOSITY)
////  void initDissipationDataStructures(struct timeStepper ts[ARRAY_ARGS 1]);
////  void destroyDissipationDataStructures(struct timeStepper ts[ARRAY_ARGS 1]);
////#endif
//
#endif /* GRIM_TIMESTEPPER_H_ */
