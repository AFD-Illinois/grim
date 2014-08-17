#include "timestepper.h"

void timeStepperInit(struct timeStepper ts[ARRAY_ARGS 1])
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
                 &ts->X1Start, &ts->X2Start, &ts->X3Start,
                 &ts->X1Size, &ts->X2Size, &ts->X3Size);
  
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVecHalfStep);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->divFluxPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->conservedVarsPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->sourceTermsPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->residualPetscVec);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->primPetscVec);

  VecSet(ts->primPetscVecOld, 0.);
  VecSet(ts->primPetscVecHalfStep, 0.);
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

  ts->timeStepCounter = 0;
  ts->tDump = 0.;

  ts->computeOldSourceTermsAndOldDivOfFluxes = 0;
  ts->computeDivOfFluxAtTimeN = 0;
  ts->computeDivOfFluxAtTimeNPlusHalf = 0;
  ts->computeSourceTermsAtTimeN = 0;
  ts->computeSourceTermsAtTimeNPlusHalf = 0;

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD, "##############################\n");
  PetscPrintf(PETSC_COMM_WORLD, "# Memory allocation complete #\n");
  PetscPrintf(PETSC_COMM_WORLD, "##############################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");

#if (RESTART)
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  if (rank==0)
  {
    if (access(RESTART_FILE, F_OK) != -1)
    {
      /* File exists */
      PetscPrintf(PETSC_COMM_WORLD, "\n");
      PetscPrintf(PETSC_COMM_WORLD, "Found restart file %s\n", RESTART_FILE); 
    }
    else
    {
      /* File does not exist */
      PetscPrintf(PETSC_COMM_WORLD, "\n");
      SETERRQ(PETSC_COMM_WORLD, 1, "Restart file %s does not exist\n",
              RESTART_FILE);
    }
  }

  PetscViewer viewer;
  PetscViewerHDF5Open(PETSC_COMM_WORLD,"restartfile.h5",
                      FILE_MODE_READ, &viewer);
  PetscObjectSetName((PetscObject) ts->primPetscVecOld, "primVars");
  VecLoad(ts->primPetscVecOld, viewer);
  PetscViewerDestroy(&viewer);

#else

  /* Set initialConditions from problem */
  initialConditions(ts);

#endif /* RESTART option */
}

/* Explicit time stepping:
 * -----------------------
 * The time stepping is done in two stages:
 *
 * 1) Get prim at t=n+1/2 using:
 *    
 *    (U^(n+1/2) - U^n)/(.5*dt) + grad(F)^(n-1) + sources^(n-1) = 0
 *    
 * 2) Now using prim at t=n+1/2, compute fluxes and sources and go from t=n to
 *    t=n+1:
 *
 *    (U^(n+1) - U^n)/dt + grad(F)^(n+1/2) + sources^(n+1/2) = 0
 *
 * IMEX time stepping:
 * -------------------
 * The time stepping is done in two stages:
 *
 * 1) Get prim at t=n+1/2 using:
 *    
 *    (U^(n+1/2) - U^n)/(.5*dt) + grad(F)^(n-1) + sources^(n-1) = 0
 *    
 * 2) Now using prim at t=n+1/2, compute fluxes but compute sources as a
 *    combination of sources at t=n and t=n+1 and go from t=n to t=n+1:
 *
 *    (U^(n+1) - U^n)/dt + grad(F)^(n+1/2) + 0.5*(sources^n + sources^(n+1)) = 0
 * 
 * Implicit time stepping:
 * -----------------------
 * Go directly from t=n to t=n+1 using
 * 
 *    (U^(n+1) - U^n)/dt + 0.5*(grad(F)^n + grad(F)^(n+1))
 *                       + 0.5*(sources^n + sources^(n+1)) = 0
 */
void timeStep(struct timeStepper ts[ARRAY_ARGS 1])
{

  PetscPrintf(PETSC_COMM_WORLD, "Step %d, time = %f\n", 
              ts->timeStepCounter, ts->t);

  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)

    /* First go from t=n to t=n+1/2 */
    ts->computeOldSourceTermsAndOldDivOfFluxes = 1;
    ts->computeDivOfFluxAtTimeN = 1;
    ts->computeDivOfFluxAtTimeNPlusHalf = 0;
    ts->computeSourceTermsAtTimeN = 1;
    ts->computeSourceTermsAtTimeNPlusHalf = 0;

    REAL dt = ts->dt;
    ts->dt = dt/2.;

    VecCopy(ts->primPetscVecOld, ts->primPetscVecHalfStep);
    SNESSolve(ts->snes, NULL, ts->primPetscVecHalfStep);

    /* Current state:
    * ts->primPetscVecHalfStep has the primitive variables at t = n+1/2 
    * ts->primPetscVecOld has the primitive variables at t = n
    */

    /* Now go from t=n to t=n+1 using
    * 1) EXPLICIT time stepping:
    * fluxes and sources computed using primitive variables at t=n+1/2.
    * 
    * 2) IMEX time stepping
    * fluxes computed using primitive variables at t=n+1/2 whereas sources
    * computed using 0.5*(Sources^n + Sources^(n-1))*/
    ts->computeOldSourceTermsAndOldDivOfFluxes = 1;
    ts->computeDivOfFluxAtTimeN = 0;
    ts->computeDivOfFluxAtTimeNPlusHalf = 1;

    #if (TIME_STEPPING==EXPLICIT)
      ts->computeSourceTermsAtTimeN = 0;
      ts->computeSourceTermsAtTimeNPlusHalf = 1;
    #elif (TIME_STEPPING==IMEX)
      ts->computeSourceTermsAtTimeN = 0;
      ts->computeSourceTermsAtTimeNPlusHalf = 0;
    #endif

    ts->dt = dt;

    VecCopy(ts->primPetscVecHalfStep, ts->primPetscVec);
    SNESSolve(ts->snes, NULL, ts->primPetscVec);

  #elif (TIME_STEPPING==IMPLICIT)

    ts->computeOldSourceTermsAndOldDivOfFluxes = 1;
    ts->computeDivOfFluxAtTimeN = 1;
    ts->computeDivOfFluxAtTimeNPlusHalf = 0;
    ts->computeSourceTermsAtTimeN = 1;
    ts->computeSourceTermsAtTimeNPlusHalf = 0;

    VecCopy(ts->primPetscVecOld, ts->primPetscVec);
    SNESSolve(ts->snes, NULL, ts->primPetscVec);

  #endif

  ts->t = ts->t + ts->dt;

  diagnostics(ts);
  problemDiagnostics(ts);

  ts->timeStepCounter++;
}

void timeStepperDestroy(struct timeStepper ts[ARRAY_ARGS 1])
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

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "# Memory deallocation complete #\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");
}
