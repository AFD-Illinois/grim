#ifndef REAPER_TIMESTEPPER_H_
#define REAPER_TIMESTEPPER_H_

#include <petsc.h>
#include "../inputs.h"
#include "../geometry/geometry.h"
#include "../gridzone/gridzone.h"
#include "../boundary/boundary.h"
#include "../reconstruct/reconstruct.h"

#define DT (.01)
#define START_TIME (0.)

#define EXPLICIT (0)
#define IMPLICIT (1)
#define IMEX     (2)

#define LOOP_OVER_TILES(X1Size, X2Size) \
for (int jTile=0; jTile<(X2Size)/TILE_SIZE_X2; jTile++) \
  for (int iTile=0; iTile<(X1Size)/TILE_SIZE_X1; iTile++)

#if (COMPUTE_DIM==2)

#define LOOP_INSIDE_TILE(iStart, iEnd, jStart, jEnd) \
for (int jInTile=jStart; jInTile<jEnd; jInTile++) \
  for (int iInTile=iStart; iInTile<iEnd; iInTile++)

#elif (COMPUTE_DIM==1)

#define LOOP_INSIDE_TILE(iStart, iEnd, jStart, jEnd) \
for (int iInTile=iStart, jInTile=0; iInTile<iEnd; iInTile++)

#endif /* LOOP_INSIDE_TILE for different COMPUTE_DIM */


struct timeStepper
{
  REAL t, dt;
  
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

/* User functions */

void timeStepperInit(struct timeStepper *ts);

void timeStep(struct timeStepper *ts);

void timeStepperDestroy(struct timeStepper *ts);

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
#endif /* REAPER_TIMESTEPPER_H_ */
