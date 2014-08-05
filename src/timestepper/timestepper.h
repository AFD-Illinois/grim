#ifndef REAPER_TIMESTEPPER_H_
#define REAPER_TIMESTEPPER_H_

#include <petsc.h>
#include "../inputs.h"
#include "../geometry/geometry.h"
#include "../gridzone/gridzone.h"
#include "../boundary/boundary.h"
#include "../reconstruct/reconstruct.h"
#include "../reaper.h"

#define EXPLICIT (0)
#define IMPLICIT (1)
#define IMEX     (2)

struct timeStepper
{
  REAL t, dt;
  
  SNES snes;
  DM dmdaWithGhostZones;
  DM dmdaWithoutGhostZones;

  Vec primPetscVec;
  Vec residualPetscVec;
  Vec primPetscVecOld;
  Vec sourceTermsPetscVecOld;
  Vec divFluxPetscVecOld;
  Vec conservedVarsPetscVecOld;

  int computedOldSourceTermsAndOldFluxes;

  int x1Start, x1Size;
  int x2Start, x2Size;
  int x3Start, x3Size;
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

#endif /* REAPER_TIMESTEPPER_H_ */
