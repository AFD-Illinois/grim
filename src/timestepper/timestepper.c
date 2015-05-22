#include "timestepper.h"

void timeStepperInit(struct timeStepper ts[ARRAY_ARGS 1])
{
  SNESCreate(PETSC_COMM_WORLD, &ts->snes);

  /* Periodic boundary conditions handled by Petsc since it is a global boundary
   * condition. Here we check for the boundary at the left edge. Obviously the
   * boundary at the right edge also must be PERIODIC if left edge is PERIODIC */
  #if (COMPUTE_DIM==1)

    #if (PHYSICAL_BOUNDARY_LEFT_EDGE==PERIODIC)
      DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, N1, DOF, NG, NULL,
                   &ts->dmdaWithGhostZones);
    #else
      DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, N1, DOF, NG, NULL,
                   &ts->dmdaWithGhostZones);
    #endif
  
  #elif (COMPUTE_DIM==2)

      #if (PHYSICAL_BOUNDARY_TOP_EDGE==PERIODIC && PHYSICAL_BOUNDARY_LEFT_EDGE==PERIODIC)
        DMDACreate2d(PETSC_COMM_WORLD, 
                     DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                     DMDA_STENCIL_BOX,
                     N1, N2,
                     PETSC_DECIDE, PETSC_DECIDE,
                     DOF, NG, PETSC_NULL, PETSC_NULL, &ts->dmdaWithGhostZones);

      #elif (PHYSICAL_BOUNDARY_LEFT_EDGE==PERIODIC)
        DMDACreate2d(PETSC_COMM_WORLD, 
                     DM_BOUNDARY_PERIODIC, DM_BOUNDARY_GHOSTED,
                     DMDA_STENCIL_BOX,
                     N1, N2,
                     PETSC_DECIDE, PETSC_DECIDE,
                     DOF, NG, PETSC_NULL, PETSC_NULL, &ts->dmdaWithGhostZones);

      #elif (PHYSICAL_BOUNDARY_TOP_EDGE==PERIODIC)
        DMDACreate2d(PETSC_COMM_WORLD, 
                     DM_BOUNDARY_GHOSTED, DM_BOUNDARY_PERIODIC,
                     DMDA_STENCIL_BOX,
                     N1, N2,
                     PETSC_DECIDE, PETSC_DECIDE,
                     DOF, NG, PETSC_NULL, PETSC_NULL, &ts->dmdaWithGhostZones);
      #else
        DMDACreate2d(PETSC_COMM_WORLD, 
                     DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                     DMDA_STENCIL_BOX,
                     N1, N2,
                     PETSC_DECIDE, PETSC_DECIDE,
                     DOF, NG, PETSC_NULL, PETSC_NULL, &ts->dmdaWithGhostZones);
      #endif

  #endif /* Choose dim and create dmdaWithGhostZones */

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
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVecLastStep);
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVecLambda);
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVecHalfStep);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->divFluxPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->conservedVarsPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->sourceTermsPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->residualPetscVec);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->OldresidualPetscVec);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->LambdaresidualPetscVec);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->LastStepresidualPetscVec);

  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
    DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->primPetscVec);
  #elif (TIME_STEPPING==IMPLICIT)
    DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVec);
  #endif

  VecSet(ts->primPetscVecOld, 0.);
  VecSet(ts->primPetscVecLastStep, 0.);
  VecSet(ts->primPetscVecLambda, 0.);
  VecSet(ts->primPetscVecHalfStep, 0.);
  VecSet(ts->divFluxPetscVecOld, 0.);
  VecSet(ts->conservedVarsPetscVecOld, 0.);
  VecSet(ts->sourceTermsPetscVecOld, 0.);
  VecSet(ts->residualPetscVec, 0.);
  VecSet(ts->OldresidualPetscVec, 0.);
  VecSet(ts->LambdaresidualPetscVec, 0.);
  VecSet(ts->LastStepresidualPetscVec, 0.);
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

  /* Now create dmda for dt */
  #if (COMPUTE_DIM==1)
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM, 0, NULL,
                 &ts->dmdaDt);
  #elif (COMPUTE_DIM==2)
    DMDACreate2d(PETSC_COMM_WORLD, 
                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX,
                 N1, N2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 COMPUTE_DIM, 0, PETSC_NULL, PETSC_NULL,
                 &ts->dmdaDt);
  #endif /* Choose dimension */

  DMCreateGlobalVector(ts->dmdaDt, &ts->dtPetscVec);

  #if (CONDUCTION)
  initConductionDataStructures(ts);
  #endif

  #if (VISCOSITY)
  initViscosityDataStructures(ts);
  #endif


  /* Initialize problem dependent data */
  PetscMalloc1(1, &ts->problemSpecificData);

  if (ts->X1Size % TILE_SIZE_X1 != 0)
  {
    PetscPrintf(PETSC_COMM_WORLD,
                "TILE_SIZE_X1 = %d does not divide X1Size = %d exactly\n",
                TILE_SIZE_X1, ts->X1Size
               );
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  #if (COMPUTE_DIM==2)
    if (ts->X2Size % TILE_SIZE_X2 != 0)
    {
      PetscPrintf(PETSC_COMM_WORLD,
                  "TILE_SIZE_X2 = %d does not divide X2Size = %d exactly\n",
                  TILE_SIZE_X2, ts->X2Size
                 );
    }
    MPI_Abort(PETSC_COMM_WORLD, 1);
  #endif

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD,
              "#################################################\n");
  PetscPrintf(PETSC_COMM_WORLD,
              "           Memory allocation complete\n\n");

  int numProcs;
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);
  PetscPrintf(PETSC_COMM_WORLD,
              " Number of MPI procs being used       = %d\n",
              numProcs);
  #if (COMPUTE_DIM==1)
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid size                            = %d\n",
                N1);
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each MPI process      = %d\n",
                ts->X1Size);
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each tile             = %d\n",
                TILE_SIZE_X1);
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
                " Grid points in each tile             = %d x %d\n",
                TILE_SIZE_X1, TILE_SIZE_X2);
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of tiles in each MPI process  = %d x %d\n",
                ts->X1Size/TILE_SIZE_X1, ts->X2Size/TILE_SIZE_X2);
  #endif

  PetscPrintf(PETSC_COMM_WORLD,
              "#################################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");

  /* Precompute the Chritoffel symbols gammaUpDownDown */
  PetscPrintf(PETSC_COMM_WORLD, "Computing Christoffel symbols...");
  setChristoffelSymbols(ts);
  PetscPrintf(PETSC_COMM_WORLD, "done\n");

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
        //SETERRQ1(PETSC_COMM_WORLD, 1, "Restart file %s does not exist\n",
        //         RESTART_FILE);
      }
    }

    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, "restartfile.h5",
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

#if (CONDUCTION)
void initConductionDataStructures(struct timeStepper ts[ARRAY_ARGS 1])
{

  #if (COMPUTE_DIM==1)
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM, 0, NULL,
                 &ts->gradTDM);

    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM*NDIM, 0, NULL,
                 &ts->graduConDM);

    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM, 0, NULL,
                 &ts->graduConHigherOrderTermsDM);

  #elif (COMPUTE_DIM==2)
    DMDACreate2d(PETSC_COMM_WORLD, 
                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX,
                 N1, N2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 COMPUTE_DIM, 0, PETSC_NULL, PETSC_NULL,
                 &ts->gradTDM);

    DMDACreate2d(PETSC_COMM_WORLD, 
                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX,
                 N1, N2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 COMPUTE_DIM*NDIM, 0, PETSC_NULL, PETSC_NULL,
                 &ts->graduConDM);

    DMDACreate2d(PETSC_COMM_WORLD, 
                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX,
                 N1, N2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 COMPUTE_DIM, 0, PETSC_NULL, PETSC_NULL,
                 &ts->graduConHigherOrderTermsDM);

  #endif

  DMCreateGlobalVector(ts->gradTDM, &ts->gradTPetscVec);
  DMCreateGlobalVector(ts->graduConDM, &ts->graduConPetscVec);
  DMCreateGlobalVector(ts->graduConHigherOrderTermsDM, 
                       &ts->graduConHigherOrderTerm1PetscVec);
  DMCreateGlobalVector(ts->graduConHigherOrderTermsDM, 
                       &ts->graduConHigherOrderTerm2PetscVec);

  VecSet(ts->gradTPetscVec, 0.);
  VecSet(ts->graduConPetscVec, 0.);
  VecSet(ts->graduConHigherOrderTerm1PetscVec, 0.);
  VecSet(ts->graduConHigherOrderTerm2PetscVec, 0.);
}

void destroyConductionDataStructures(struct timeStepper ts[ARRAY_ARGS 1])
{
  VecDestroy(&ts->gradTPetscVec);
  VecDestroy(&ts->graduConPetscVec);
  VecDestroy(&ts->graduConHigherOrderTerm1PetscVec);
  VecDestroy(&ts->graduConHigherOrderTerm2PetscVec);

  DMDestroy(&ts->gradTDM);
  DMDestroy(&ts->graduConDM);
  DMDestroy(&ts->graduConHigherOrderTermsDM);
}
#endif


#if (VISCOSITY)
void initViscosityDataStructures(struct timeStepper ts[ARRAY_ARGS 1])
{

  #if (COMPUTE_DIM==1)
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM*NDIM, 0, NULL,
                 &ts->graduConVisDM);
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM, 0, NULL,
                 &ts->graduConHigherOrderTermsVisDM);

  #elif (COMPUTE_DIM==2)
    DMDACreate2d(PETSC_COMM_WORLD, 
                   DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                   DMDA_STENCIL_BOX,
                   N1, N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   COMPUTE_DIM*NDIM, 0, PETSC_NULL, PETSC_NULL,
                   &ts->graduConVisDM);
    DMDACreate2d(PETSC_COMM_WORLD, 
                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX,
                 N1, N2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 COMPUTE_DIM, 0, PETSC_NULL, PETSC_NULL,
                 &ts->graduConHigherOrderTermsVisDM);

  #endif

  DMCreateGlobalVector(ts->graduConVisDM, &ts->graduConVisPetscVec);
  DMCreateGlobalVector(ts->graduConHigherOrderTermsVisDM, 
                       &ts->graduConHigherOrderTerm1VisPetscVec);
  DMCreateGlobalVector(ts->graduConHigherOrderTermsVisDM, 
                       &ts->graduConHigherOrderTerm2VisPetscVec);

  VecSet(ts->graduConVisPetscVec, 0.);
  VecSet(ts->graduConHigherOrderTerm1VisPetscVec, 0.);
  VecSet(ts->graduConHigherOrderTerm2VisPetscVec, 0.);
}

void destroyViscosityDataStructures(struct timeStepper ts[ARRAY_ARGS 1])
{
  VecDestroy(&ts->graduConVisPetscVec);
  VecDestroy(&ts->graduConHigherOrderTerm1VisPetscVec);
  VecDestroy(&ts->graduConHigherOrderTerm2VisPetscVec);

  DMDestroy(&ts->graduConVisDM);
  DMDestroy(&ts->graduConHigherOrderTermsVisDM);
}
#endif

void setChristoffelSymbols(struct timeStepper ts[ARRAY_ARGS 1])
{
  VecSet(ts->connectionPetscVec, 0.);

  ARRAY(connectionGlobal);
  DMDAVecGetArrayDOF(ts->connectionDMDA, ts->connectionPetscVec,
                     &connectionGlobal);

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
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
  PetscPrintf(PETSC_COMM_WORLD, "Start time step \n");
  PetscErrorCode errCode;

  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)

    /* First go from t=n to t=n+1/2 */
    ts->computeOldSourceTermsAndOldDivOfFluxes = 1;
    ts->computeDivOfFluxAtTimeN = 1;
    ts->computeDivOfFluxAtTimeNPlusHalf = 0;
    ts->computeSourceTermsAtTimeN = 1;
    ts->computeSourceTermsAtTimeNPlusHalf = 0;

    REAL dt = ts->dt;
    ts->dt = dt/2.;

    int MAX_IT = 3;
    int MAX_LAMBDA_LOOPS = 10;
    double ALPHA_LS = 1.e-4;
    double LAMBDA_MIN = 1.e-12;
    REAL atol = 1.e-7;
    REAL rtol = 1.e-5;

    ARRAY(OldresidualLocal);
    ARRAY(LambdaresidualLocal);
    ARRAY(LastStepresidualLocal);
    ARRAY(residualLocal);
    ARRAY(HalfStepprimVecLocal);
    ARRAY(OldprimVecLocal);
    ARRAY(LastStepprimVecLocal);
    ARRAY(primVecLocal);
    ARRAY(LambdaprimVecLocal);
    
    
    int Npoints = ts->X1Size*ts->X2Size;
    int UseCustomLineSearch = 1;
    REAL mLambda[Npoints];
    int mConverged[Npoints];

    if(UseCustomLineSearch==1)
      {
	errCode = computeResidual(ts->snes,ts->primPetscVecOld,ts->OldresidualPetscVec,ts);
	VecCopy(ts->primPetscVecOld, ts->primPetscVecHalfStep);
	VecCopy(ts->OldresidualPetscVec, ts->residualPetscVec);
	DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
			   &OldresidualLocal);
	DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecOld,
			   &OldprimVecLocal);
	//Newton linesearch method
	for(int it=0;it<MAX_IT;it++)
	  {
	    VecCopy(ts->residualPetscVec,ts->LastStepresidualPetscVec);
	    VecCopy(ts->primPetscVecHalfStep, ts->primPetscVecLastStep);
	    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecLastStep,
			       &LastStepprimVecLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LastStepresidualPetscVec,
			       &LastStepresidualLocal);
	    //1) Standard petsc inversion. Should use basic linesearch, and snes_max_it=1
	    errCode = SNESSolve(ts->snes, NULL, ts->primPetscVecHalfStep);
	    int AllPointsConverged = 0;
	    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
			       &primVecLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
			       &residualLocal);
	    VecCopy(ts->residualPetscVec,ts->LambdaresidualPetscVec);
	    VecCopy(ts->primPetscVecHalfStep, ts->primPetscVecLambda);
	    //Linesearch within each step
	    for(int itlam = 0; itlam<MAX_LAMBDA_LOOPS && AllPointsConverged==0;itlam++)
	      {
		DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecLambda,
				   &LambdaprimVecLocal);
		DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LambdaresidualPetscVec,
				   &LambdaresidualLocal);
		//Compute L2 norm of residual pointwise
		AllPointsConverged = 1;
		int NptsNotConv = 0;
		REAL DiagOldRes = 0.;
		REAL DiagNewRes = 0.;
                #if (USE_OPENMP)
                  #pragma omp parallel for
                #endif
		LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
		  {
		    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
		      {
			int idx = iTile*ts->X2Size*TILE_SIZE_X1
			  +jTile*TILE_SIZE_X1*TILE_SIZE_X2+
			  iInTile*TILE_SIZE_X2+jInTile;
			if(itlam==0)
			  {
			    mLambda[idx]=1.;
			    mConverged[idx]=0;
			  }
			struct gridZone zone;
			setGridZone(iTile, jTile,
				    iInTile, jInTile,
				    ts->X1Start, ts->X2Start,
				    ts->X1Size, ts->X2Size,
				    &zone);
			REAL Oldres = 0.;
			REAL Lambdares = 0.;
			REAL Firstres = 0.;
			REAL temp = 0.;
			for (int var=0; var<DOF; var++)
			  {
			    temp = INDEX_PETSC(LastStepresidualLocal,&zone,var);
			    Oldres+=temp*temp;
			    temp = INDEX_PETSC(LambdaresidualLocal,&zone,var);
			    Lambdares+=temp*temp; 
			    temp = INDEX_PETSC(OldresidualLocal,&zone,var);
			    Firstres+=temp*temp;
			  }
			DiagOldRes += Oldres;
			DiagNewRes += Lambdares;
			//Stop if we already converged (we went this far just
			//to include all points in the residual)
			if(mConverged[idx]==1)
			  continue;
			if(Lambdares<Oldres*(1.-mLambda[idx]*ALPHA_LS))
			  {
			    //Acceptable end to the linesearch
			    mConverged[idx]=1;
			  }
			else
			  {
			    //if the residual was already small, just accept it...
			    if(sqrt(Oldres)<atol+rtol*sqrt(Firstres))
			      {
				for (int var=0; var<DOF; var++)
				  INDEX_PETSC(LambdaprimVecLocal,&zone,var)= 
				    INDEX_PETSC(LastStepprimVecLocal,&zone,var);
				mConverged[idx]=1;
				mLambda[idx]=0.;
				DiagNewRes+=Oldres-Lambdares;
			      }
			    else
			      {
				//backtrack
				REAL lam = mLambda[idx];
				REAL lamupd = Oldres*lam*lam/(Lambdares+(2.*lam-1.)*Oldres);
				if(lamupd>0.5*lam)
				  lamupd=0.5*lam;
				if(lamupd<0.1*lam)
				  lamupd=0.1*lam;
				mLambda[idx]=lamupd;
				if(mLambda[idx]<LAMBDA_MIN)
				  mLambda[idx]=LAMBDA_MIN;
				for (int var=0; var<DOF; var++)
				  INDEX_PETSC(LambdaprimVecLocal,&zone,var)=
				    (1.-mLambda[idx])*INDEX_PETSC(LastStepprimVecLocal,&zone,var)+
				    mLambda[idx]*INDEX_PETSC(primVecLocal,&zone,var);
				AllPointsConverged=0;
				NptsNotConv++;
			      }
			  }
		      }
		  }
		DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
				       ts->primPetscVecLambda,
				       &LambdaprimVecLocal);
		DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, 
				       ts->LambdaresidualPetscVec,
				       &LambdaresidualLocal);
		DiagOldRes = sqrt(DiagOldRes);
		DiagNewRes = sqrt(DiagNewRes);
		if(itlam==0)
		  printf("Step %i - Initial residual = %e \n",it,DiagOldRes);
		printf("New residual = %e. ",DiagNewRes);
		printf("%i of %i points have not converged.\n",NptsNotConv,Npoints);
		errCode = computeResidual(ts->snes,ts->primPetscVecLambda,ts->LambdaresidualPetscVec,ts);
	      }
	    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
				   ts->primPetscVecHalfStep,
				   &primVecLocal);
	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, 
				   ts->residualPetscVec,
				   &residualLocal);
	    VecCopy(ts->LambdaresidualPetscVec,ts->residualPetscVec);
	    VecCopy(ts->primPetscVecLambda, ts->primPetscVecHalfStep);
	    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
				   ts->primPetscVecLastStep,
				   &LastStepprimVecLocal);
	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, 
				   ts->LastStepresidualPetscVec,
				   &LastStepresidualLocal);
	  }
	DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
			       &OldresidualLocal);
	DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecOld,
			       &OldprimVecLocal);
      }
    else
      {
	VecCopy(ts->primPetscVecOld, ts->primPetscVecHalfStep);
	errCode = SNESSolve(ts->snes, NULL, ts->primPetscVecHalfStep);
      }
    //CHKERRQ(errCode);

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

    if(UseCustomLineSearch==1)
      {
	errCode = computeResidual(ts->snes,ts->primPetscVecHalfStep,ts->OldresidualPetscVec,ts);
	VecCopy(ts->primPetscVecHalfStep, ts->primPetscVec);
	VecCopy(ts->OldresidualPetscVec, ts->residualPetscVec);
	DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
			   &OldresidualLocal);
	DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
			   &HalfStepprimVecLocal);
	for(int it=0;it<MAX_IT;it++)
	  {
	    VecCopy(ts->residualPetscVec,ts->LastStepresidualPetscVec);
	    VecCopy(ts->primPetscVec, ts->primPetscVecLastStep);
	    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecLastStep,
			       &LastStepprimVecLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LastStepresidualPetscVec,
			       &LastStepresidualLocal);
	    errCode = SNESSolve(ts->snes, NULL, ts->primPetscVec);
	    int AllPointsConverged = 0;
	    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVec,
			       &primVecLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
			       &residualLocal);
	    VecCopy(ts->residualPetscVec,ts->LambdaresidualPetscVec);
	    VecCopy(ts->primPetscVec, ts->primPetscVecLambda);
	    //Linesearch within each step
	    for(int itlam = 0; itlam<MAX_LAMBDA_LOOPS && AllPointsConverged==0;itlam++)
	      {
		DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecLambda,
				   &LambdaprimVecLocal);
		DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LambdaresidualPetscVec,
				   &LambdaresidualLocal);
		//Compute L2 norm of residual pointwise                                                                                                                         
		AllPointsConverged = 1;
		int NptsNotConv = 0;
		REAL DiagOldRes = 0.;
		REAL DiagNewRes = 0.;
                #if (USE_OPENMP)
                  #pragma omp parallel for
                #endif
		LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
		  {
		    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
		      {
			int idx = iTile*ts->X2Size*TILE_SIZE_X1
			  +jTile*TILE_SIZE_X1*TILE_SIZE_X2+
			  iInTile*TILE_SIZE_X2+jInTile;
			if(itlam==0)
			  {
			    mLambda[idx]=1.;
			    mConverged[idx]=0;
			  }
			struct gridZone zone;
			setGridZone(iTile, jTile,
				    iInTile, jInTile,
				    ts->X1Start, ts->X2Start,
				    ts->X1Size, ts->X2Size,
				    &zone);
			REAL Oldres = 0.;
			REAL Lambdares = 0.;
			REAL Firstres = 0.;
			REAL temp = 0.;
			for (int var=0; var<DOF; var++)
			  {
			    temp = INDEX_PETSC(OldresidualLocal,&zone,var);
			    Firstres+=temp*temp;
			    temp = INDEX_PETSC(LastStepresidualLocal,&zone,var);
			    Oldres+=temp*temp;
			    temp = INDEX_PETSC(LambdaresidualLocal,&zone,var);
			    Lambdares+=temp*temp; 
			  }
			DiagOldRes += Oldres;
			DiagNewRes += Lambdares;
			//Stop if we already converged (we went this far just
			//to include all points in the residual)
			if(mConverged[idx]==1)
			  continue;
			if(Lambdares<Oldres*(1.-mLambda[idx]*ALPHA_LS))
			  {
			    //Acceptable end to the linesearch
			    mConverged[idx]=1;
			  }
			else
			  {
			    //if the residual was already small, just accept it...
			    if(sqrt(Oldres)<atol+rtol*sqrt(Firstres))
			      {
				for (int var=0; var<DOF; var++)
				  INDEX_PETSC(LambdaprimVecLocal,&zone,var)= 
				    INDEX_PETSC(LastStepprimVecLocal,&zone,var);
				mConverged[idx]=1;
				mLambda[idx]=0.;
				DiagNewRes+=Oldres-Lambdares;
			      }
			    else
			      {
				//backtrack
				REAL lam = mLambda[idx];
				REAL lamupd = Oldres*lam*lam/(Lambdares+(2.*lam-1.)*Oldres);
				if(lamupd>0.5*lam)
                                  lamupd=0.5*lam;
                                if(lamupd<0.1*lam)
                                  lamupd=0.1*lam;
                                mLambda[idx]=lamupd;
				if(mLambda[idx]<LAMBDA_MIN)
				  mLambda[idx]=LAMBDA_MIN;
				for (int var=0; var<DOF; var++)
				  INDEX_PETSC(LambdaprimVecLocal,&zone,var)=
				    (1.-mLambda[idx])*INDEX_PETSC(LastStepprimVecLocal,&zone,var)+
				    mLambda[idx]*INDEX_PETSC(primVecLocal,&zone,var);
				AllPointsConverged=0;
				NptsNotConv++;
			      }
			  }
		      }
		  }
		DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
				       ts->primPetscVecLambda,
				       &LambdaprimVecLocal);
		DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
				       ts->LambdaresidualPetscVec,
				       &LambdaresidualLocal);
		DiagOldRes = sqrt(DiagOldRes);
		DiagNewRes = sqrt(DiagNewRes);
		if(itlam==0 && it==0)
		  printf("Initial residual = %e\n",DiagOldRes);
		if(itlam==0)
		  printf("Step %i \n",it);
		printf("New residual = %e. ",DiagNewRes);
		printf("%i of %i points have not converged.\n",NptsNotConv,Npoints);
		errCode = computeResidual(ts->snes,ts->primPetscVecLambda,ts->LambdaresidualPetscVec,ts);
	      }
	    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
				   ts->primPetscVec,
				   &primVecLocal);
	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, 
				   ts->residualPetscVec,
				   &residualLocal);
	    VecCopy(ts->LambdaresidualPetscVec,ts->residualPetscVec);
	    VecCopy(ts->primPetscVecLambda, ts->primPetscVec);
	    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
				   ts->primPetscVecLastStep,
				   &LastStepprimVecLocal);
	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
				   ts->LastStepresidualPetscVec,
				   &LastStepresidualLocal);
	  }
	DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
			       &OldresidualLocal);
	DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
			       &HalfStepprimVecLocal);
      }    
    else
      {
	VecCopy(ts->primPetscVecHalfStep, ts->primPetscVec);
	errCode = SNESSolve(ts->snes, NULL, ts->primPetscVec);
      }
    //Diagnostics
    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVec,
		       &primVecLocal);
    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
		       &residualLocal);   
    int NptsResMag[16];
    for(int i=0;i<16;i++)
      NptsResMag[i]=0;
    #if (USE_OPENMP)
      #pragma omp parallel for
    #endif
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
	    REAL res = 0.;
	    REAL temp = 0.;
	    for (int var=0; var<DOF; var++)
	      {
		temp = INDEX_PETSC(residualLocal,&zone,var);
		res+=temp*temp;
	      }
	    res=sqrt(res);
	    res=log10(res+1.e-15);
	    for(int i=0;i<16;i++)
	      if(res>-i)
		{
		  NptsResMag[i]++;
		  break;
		}
	  }
      }
    printf("Distribution of residuals : ");
    for(int i=0;i<16;i++)
      printf("%i ;", NptsResMag[i]);
    printf("\n");
    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
			   ts->primPetscVec,
			   &primVecLocal);
    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
			   ts->residualPetscVec,
			   &residualLocal);    
    //CHKERRQ(errCode);

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
  ts->timeStepCounter++;
  PetscPrintf(PETSC_COMM_WORLD,
              "\nCompleted step %d, current time = %.5f, dt = %.5f\n\n",
              ts->timeStepCounter, ts->t, ts->dt);

  /* Problem dependent full step diagnostics */
  fullStepDiagnostics(ts);

  /* Copy solved variables to old variables */
  VecCopy(ts->primPetscVec, ts->primPetscVecOld);

  /* Finally, take a dump when needed */
  diagnostics(ts);

  REAL newDt;
  #if (COMPUTE_DIM==1)

    VecStrideMin(ts->dtPetscVec, 0, NULL, &newDt);

  #elif (COMPUTE_DIM==2)
    REAL dtX1, dtX2;
    VecStrideMin(ts->dtPetscVec, 0, NULL, &dtX1);
    VecStrideMin(ts->dtPetscVec, 1, NULL, &dtX2);

    newDt = 1./( (1./dtX1) + (1./dtX2) );
  #endif

  if (newDt > MAX_DT_INCREMENT*ts->dt)
  {
    newDt = MAX_DT_INCREMENT*ts->dt;
  }

  ts->dt = newDt;
}

void timeStepperDestroy(struct timeStepper ts[ARRAY_ARGS 1])
{
  VecDestroy(&ts->primPetscVecOld);
  VecDestroy(&ts->primPetscVecLastStep);
  VecDestroy(&ts->primPetscVecLambda);
  VecDestroy(&ts->divFluxPetscVecOld);
  VecDestroy(&ts->conservedVarsPetscVecOld);
  VecDestroy(&ts->sourceTermsPetscVecOld);
  VecDestroy(&ts->residualPetscVec);
  VecDestroy(&ts->OldresidualPetscVec);
  VecDestroy(&ts->LastStepresidualPetscVec);
  VecDestroy(&ts->LambdaresidualPetscVec);
  VecDestroy(&ts->primPetscVec);
  VecDestroy(&ts->connectionPetscVec);

  DMDestroy(&ts->dmdaWithGhostZones);
  DMDestroy(&ts->dmdaWithoutGhostZones);
  DMDestroy(&ts->connectionDMDA);

  SNESDestroy(&ts->snes);

  PetscFree(ts->problemSpecificData);

  #if (CONDUCTION)
    destroyConductionDataStructures(ts);
  #endif

  #if (VISCOSITY)
    destroyViscosityDataStructures(ts);
  #endif

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "# Memory deallocation complete #\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");
}
