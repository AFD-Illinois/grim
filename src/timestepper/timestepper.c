#include "timestepper.h"

void timeStepperInit(struct timeStepper *ts)
{
  SNESCreate(PETSC_COMM_WORLD, &ts->snes);

/* Periodic boundary conditions handled by Petsc since it is a global boundary
 * condition. Here we check for the boundary at the left edge. Obviously the
 * boundary at the right edge also must be PERIODIC if left edge is PERIODIC */
#if (PHYSICAL_BOUNDARY_LEFT_EDGE==PERIODIC)
#if (COMPUTE_DIM==1)
  DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, N1, DOF, NG, NULL,
               &ts->dmdaWithGhostZones);
#elif (COMPUTE_DIM==2)
  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
               DMDA_STENCIL_BOX,
               N1, N2,
               PETSC_DECIDE, PETSC_DECIDE,
               DOF, NG, PETSC_NULL, PETSC_NULL, &ts->dmdaWithGhostZones);
#endif /* Choose dimension */
#else
#if (COMPUTE_DIM==1)
  DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, N1, DOF, NG, NULL,
               &ts->dmdaWithGhostZones);
#elif (COMPUTE_DIM==2)
  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
               DMDA_STENCIL_BOX,
               N1, N2,
               PETSC_DECIDE, PETSC_DECIDE,
               DOF, NG, PETSC_NULL, PETSC_NULL, &ts->dmdaWithGhostZones);
#endif /* Choose dimension */
#endif /* Create mesh structure */

/* If explicit or imex time stepping, then create another DM with no ghost
 * zones. Graph coloring with ensure that the number of evaluations to construct
 * the jacobian will be equal to DOF^2 */
#if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
#if (COMPUTE_DIM==1)
  DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, DOF, 0, NULL,
               &ts->dmdaWithoutGhostZones);
#elif (COMPUTE_DIM==2)
  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
               DMDA_STENCIL_BOX,
               N1, N2,
               PETSC_DECIDE, PETSC_DECIDE,
               DOF, 0, PETSC_NULL, PETSC_NULL, &ts->dmdaWithoutGhostZones);
#endif /* Choose dimension */
  SNESSetDM(ts->snes, ts->dmdaWithoutGhostZones);
#elif (TIME_STEPPING==IMPLICIT)
  /* If implicit time stepping, then use the DM created earlier which has ghost
   * zones, in SNES */
  SNESSetDM(ts->snes, ts->dmdaWithGhostZones);
#endif

  DMDAGetCorners(ts->dmdaWithGhostZones,
                 &ts->x1Start, &ts->x2Start, &ts->x3Start,
                 &ts->x1Size, &ts->x2Size, &ts->x3Size);
  
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->divFluxPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->conservedVarsPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->sourceTermsPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->residualPetscVec);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->primPetscVec);

  VecSet(ts->primPetscVecOld, 0.);
  VecSet(ts->divFluxPetscVecOld, 0.);
  VecSet(ts->conservedVarsPetscVecOld, 0.);
  VecSet(ts->sourceTermsPetscVecOld, 0.);
  VecSet(ts->residualPetscVec, 0.);
  VecSet(ts->primPetscVec, 0.);

  SNESSetFunction(ts->snes, ts->residualPetscVec,
                  computeResidual, ts);

  SNESSetFromOptions(ts->snes);

  ts->dt = DT;
  ts->t = START_TIME;

  ts->computedOldSourceTermsAndOldFluxes = 0;

  PetscPrintf(PETSC_COMM_WORLD, "|--------------------------|\n");
  PetscPrintf(PETSC_COMM_WORLD, "|Memory allocation complete|\n");
  PetscPrintf(PETSC_COMM_WORLD, "|--------------------------|\n");
}

void timeStep(struct timeStepper *ts)
{
  ts->computedOldSourceTermsAndOldFluxes = 0;

  VecCopy(ts->primPetscVecOld, ts->primPetscVec);
  SNESSolve(ts->snes, NULL, ts->primPetscVec);

  ts->t = ts->t + ts->dt;
}

void timeStepperDestroy(struct timeStepper *ts)
{
  VecDestroy(&ts->primPetscVecOld);
  VecDestroy(&ts->divFluxPetscVecOld);
  VecDestroy(&ts->conservedVarsPetscVecOld);
  VecDestroy(&ts->sourceTermsPetscVecOld);
  VecDestroy(&ts->residualPetscVec);
  VecDestroy(&ts->primPetscVec);

  DMDestroy(&ts->dmdaWithGhostZones);
  DMDestroy(&ts->dmdaWithoutGhostZones);

  SNESDestroy(&ts->snes);

  PetscPrintf(PETSC_COMM_WORLD, "|----------------------------|\n");
  PetscPrintf(PETSC_COMM_WORLD, "|Memory deallocation complete|\n");
  PetscPrintf(PETSC_COMM_WORLD, "|----------------------------|\n");
}
