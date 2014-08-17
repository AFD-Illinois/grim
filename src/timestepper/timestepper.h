#ifndef GRIM_TIMESTEPPER_H_
#define GRIM_TIMESTEPPER_H_

#include <petsc.h>
#include <petscviewerhdf5.h>
#include "../inputs.h"
#include "../geometry/geometry.h"
#include "../gridzone/gridzone.h"
#include "../boundary/boundary.h"
#include "../reconstruct/reconstruct.h"
#include "../riemannsolver/riemannsolver.h"

#define EXPLICIT (0)
#define IMPLICIT (1)
#define IMEX     (2)

struct timeStepper
{
  REAL t, dt, tDump;
  int timeStepCounter;

  SNES snes;
  DM dmdaWithGhostZones;
  DM dmdaWithoutGhostZones;

  Vec primPetscVec;
  Vec residualPetscVec;
  Vec primPetscVecOld;
  Vec primPetscVecHalfStep;
  Vec sourceTermsPetscVecOld;
  Vec divFluxPetscVecOld;
  Vec conservedVarsPetscVecOld;

  int computeOldSourceTermsAndOldDivOfFluxes;
  int computeDivOfFluxAtTimeN;
  int computeDivOfFluxAtTimeNPlusHalf;
  int computeSourceTermsAtTimeN;
  int computeSourceTermsAtTimeNPlusHalf;

  int X1Start, X1Size;
  int X2Start, X2Size;
  int X3Start, X3Size;
};

/* The problem library requires timeStepper struct, so put this include after
 * the definition of the struct */
#include "../problem/problem.h"

/* User functions */

void timeStepperInit(struct timeStepper ts[ARRAY_ARGS 1]);

void timeStep(struct timeStepper ts[ARRAY_ARGS 1]);

void timeStepperDestroy(struct timeStepper ts[ARRAY_ARGS 1]);

/* Internal functions */

PetscErrorCode computeResidual(SNES snes, 
                               Vec primPetscVec,
                               Vec residalPetscVec,
                               void *ptr);

void computeFluxesOverTile(const REAL primTile[ARRAY_ARGS TILE_SIZE],
                           const int iTile, const int jTile,
                           const int X1Start, const int X2Start,
                           const int X1Size, const int X2Size,
                           REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
                           REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]);

void diagnostics(struct timeStepper ts[ARRAY_ARGS 1]);

#endif /* GRIM_TIMESTEPPER_H_ */
