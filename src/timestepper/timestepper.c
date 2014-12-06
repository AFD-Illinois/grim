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

  #else /* Not a periodic boundary */
    
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

  #endif /* Create dmdaWithGhostZones */

  /* Now create dmdaWithoutGhostZones for the vectors that don't need
   * communication */
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

  /* If explicit or imex time stepping, then create another DM with no ghost
   * zones. Graph coloring with ensure that the number of evaluations to construct
   * the jacobian will be equal to DOF^2 */
  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
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

  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
    DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->primPetscVec);
  #elif (TIME_STEPPING==IMPLICIT)
    DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVec);
  #endif

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
  ts->tDump = START_TIME;

  ts->timeStepCounter = 0;
  ts->dumpCounter = 0;

  ts->computeOldSourceTermsAndOldDivOfFluxes = 0;
  ts->computeDivOfFluxAtTimeN = 0;
  ts->computeDivOfFluxAtTimeNPlusHalf = 0;
  ts->computeSourceTermsAtTimeN = 0;
  ts->computeSourceTermsAtTimeNPlusHalf = 0;

  /* Now create dmda for the connection coefficients */
  #if (COMPUTE_DIM==1)
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, 64, 0, NULL,
                 &ts->connectionDMDA);
  #elif (COMPUTE_DIM==2)
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX, N1, N2, PETSC_DECIDE, PETSC_DECIDE,
                 64, 0, PETSC_NULL, PETSC_NULL, &ts->connectionDMDA);
  #endif /* Choose dimension */

  DMCreateGlobalVector(ts->connectionDMDA, &ts->connectionPetscVec);

  if (ts->X1Size % TILE_SIZE_X1 != 0)
  {
    SETERRQ2(PETSC_COMM_WORLD, 1,
             "TILE_SIZE_X1 = %d does not divide X1Size = %d exactly\n",
             TILE_SIZE_X1, ts->X1Size);
  }
  #if (COMPUTE_DIM==2)
    if (ts->X2Size % TILE_SIZE_X2 != 0)
    {
      SETERRQ2(PETSC_COMM_WORLD, 1,
               "TILE_SIZE_X2 = %d does not divide X2Size = %d exactly\n",
               TILE_SIZE_X2, ts->X2Size);
    }
  #endif

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD,
              "#################################################\n");
  PetscPrintf(PETSC_COMM_WORLD,
              "           Memory allocation complete\n\n");

  int numProcs;
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);
  PetscPrintf(PETSC_COMM_WORLD,
              " Number of processors being used      = %d\n",
              numProcs);
  #if (COMPUTE_DIM==1)
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid size                            = %d\n",
                N1);
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each MPI process      = %d\n",
                ts->X1Size);
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of tiles in each MPI process  = %d\n",
                ts->X1Size/TILE_SIZE_X1);
  #elif (COMPUTE_DIM==2)
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid size                            = %d x %d\n",
                N1, N2);
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each MPI process      = %d x %d\n",
                ts->X1Size, ts->X2Size);
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of tiles in each MPI process  = %d x %d\n",
                ts->X1Size/TILE_SIZE_X1, ts->X2Size/TILE_SIZE_X2);
  #endif

  PetscPrintf(PETSC_COMM_WORLD,
              "#################################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");

  /* Precompute the Chritoffel symbols gammaUpDownDown */
  setChristoffelSymbols(ts);

  #if (RESTART)
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if (rank==0)
    {
      if (access(RESTART_FILE, F_OK) != -1)
      {
        /* File exists */
        PetscPrintf(PETSC_COMM_WORLD, "\nFound restart file: %s\n\n", RESTART_FILE); 
      }
      else
      {
        /* File does not exist */
        PetscPrintf(PETSC_COMM_WORLD, "\n");
        SETERRQ1(PETSC_COMM_WORLD, 1, "Restart file %s does not exist\n",
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

  /* Output the initial conditions */
  VecCopy(ts->primPetscVecOld, ts->primPetscVec);
  diagnostics(ts);

}

void setChristoffelSymbols(struct timeStepper ts[ARRAY_ARGS 1])
{
  VecSet(ts->connectionPetscVec, 0.);

  ARRAY(connectionGlobal);
  DMDAVecGetArrayDOF(ts->connectionDMDA, ts->connectionPetscVec,
                     &connectionGlobal);

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
      struct geometry geom; setGeometry(XCoords, &geom);

      /* Now compute connectionGlobal with 
       * Index Up   - eta
       * Index down - mu
       * Index down - nu */
      for (int eta=0; eta<NDIM; eta++)
      {
        for (int mu=0; mu<NDIM; mu++)
        {
          for (int nu=0; nu<NDIM; nu++)
          {
            for (int alpha=0; alpha<NDIM; alpha++)
            {
              #if (COMPUTE_DIM==1)
                connectionGlobal[zone.i][GAMMA_UP_DOWN_DOWN(eta, mu, nu)]
              #elif (COMPUTE_DIM==2)
                connectionGlobal[zone.j][zone.i][GAMMA_UP_DOWN_DOWN(eta, mu, nu)]
              #endif
              +=
                geom.gCon[eta][alpha]
              * gammaDownDownDown(alpha, mu, nu, XCoords);
            }
          }
        }
      }

    }
  }

  DMDAVecRestoreArrayDOF(ts->connectionDMDA, ts->connectionPetscVec,
                         &connectionGlobal);
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

    /* Problem dependent half step diagnostics */
    halfStepDiagnostics(ts);

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
    VecCopy(ts->primPetscVec, ts->primPetscVecOld);

  #elif (TIME_STEPPING==IMPLICIT)

    ts->computeOldSourceTermsAndOldDivOfFluxes = 1;
    ts->computeDivOfFluxAtTimeN = 1;
    ts->computeDivOfFluxAtTimeNPlusHalf = 0;
    ts->computeSourceTermsAtTimeN = 1;
    ts->computeSourceTermsAtTimeNPlusHalf = 0;

    VecCopy(ts->primPetscVecOld, ts->primPetscVec);
    SNESSolve(ts->snes, NULL, ts->primPetscVec);
    VecCopy(ts->primPetscVec, ts->primPetscVecOld);

  #endif

  ts->t = ts->t + ts->dt;
  ts->timeStepCounter++;

  PetscPrintf(PETSC_COMM_WORLD,
              "\nCompleted step %d, current time = %.5f\n\n",
              ts->timeStepCounter, ts->t, ts->tDump);

  /* Problem dependent full step diagnostics */
  fullStepDiagnostics(ts);

  diagnostics(ts);
}

void timeStepperDestroy(struct timeStepper ts[ARRAY_ARGS 1])
{
  VecDestroy(&ts->primPetscVecOld);
  VecDestroy(&ts->divFluxPetscVecOld);
  VecDestroy(&ts->conservedVarsPetscVecOld);
  VecDestroy(&ts->sourceTermsPetscVecOld);
  VecDestroy(&ts->residualPetscVec);
  VecDestroy(&ts->primPetscVec);
  VecDestroy(&ts->connectionPetscVec);

  DMDestroy(&ts->dmdaWithGhostZones);
  DMDestroy(&ts->dmdaWithoutGhostZones);
  DMDestroy(&ts->connectionDMDA);

  SNESDestroy(&ts->snes);

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "# Memory deallocation complete #\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");
}
