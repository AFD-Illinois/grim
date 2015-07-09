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
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVecTmp);
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVecLambda);
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVecHalfStep);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->divFluxPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->conservedVarsPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->sourceTermsPetscVecOld);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->residualPetscVec);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->OldresidualPetscVec);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->LambdaresidualPetscVec);

  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
    DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &ts->primPetscVec);
  #elif (TIME_STEPPING==IMPLICIT)
    DMCreateGlobalVector(ts->dmdaWithGhostZones, &ts->primPetscVec);
  #endif

  VecSet(ts->primPetscVecOld, 0.);
  VecSet(ts->primPetscVecTmp, 0.);
  VecSet(ts->primPetscVecLambda, 0.);
  VecSet(ts->primPetscVecHalfStep, 0.);
  VecSet(ts->divFluxPetscVecOld, 0.);
  VecSet(ts->conservedVarsPetscVecOld, 0.);
  VecSet(ts->sourceTermsPetscVecOld, 0.);
  VecSet(ts->residualPetscVec, 0.);
  VecSet(ts->OldresidualPetscVec, 0.);
  VecSet(ts->LambdaresidualPetscVec, 0.);
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
    REAL atol = 1.e-8;
    REAL rtol = 1.e-6;
    REAL OrigRes = 0.;
    ARRAY(OldresidualLocal);
    ARRAY(LambdaresidualLocal);
    ARRAY(residualLocal);
    ARRAY(OldprimVecLocal);
    ARRAY(TmpprimVecLocal);
    ARRAY(primVecLocal);
    ARRAY(LambdaprimVecLocal);
    VecCopy(ts->primPetscVecOld, ts->primPetscVecHalfStep);

    for(int it=0;it<MAX_IT;it++)
      {
	VecCopy(ts->primPetscVecHalfStep, ts->primPetscVecTmp);
	errCode = computeResidual(ts->snes,ts->primPetscVecTmp,ts->OldresidualPetscVec,ts);
	errCode = SNESSolve(ts->snes, NULL, ts->primPetscVecHalfStep);
	int AllPointsConverged = 0;
	REAL minlambda = 1.e-5;
	
	int MAX_LAMBDA_LOOPS = 3;
	DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
			   &OldresidualLocal);
	DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->primPetscVecTmp,
			   &TmpprimVecLocal);
	for(int itlam = 0; itlam<MAX_LAMBDA_LOOPS && AllPointsConverged==0;itlam++)
	  {
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
			       &residualLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LambdaresidualPetscVec,
                               &LambdaresidualLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->primPetscVecHalfStep,
			       &primVecLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->primPetscVecLambda,
                               &LambdaprimVecLocal);
	    REAL DiagOldRes = 0.;
	    REAL DiagFullRes = 0.;
	    REAL DiagLamRes = 0.;

	    AllPointsConverged = 1;
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
		    REAL Oldres = 0.;
		    REAL temp = 0.;
		    for (int var=0; var<DOF; var++)
		      {
			temp = INDEX_PETSC(OldresidualLocal,&zone,var);
			Oldres+=temp*temp;
			temp = INDEX_PETSC(residualLocal,&zone,var);
			res+=temp*temp;
		      }
		    res=sqrt(res);
		    Oldres=sqrt(Oldres);
		    DiagOldRes += Oldres;
		    DiagFullRes += res;
		    if(fabs(res)<fabs(Oldres))
		      AllPointsConverged=0;
		    REAL lambda = 1.;
		    if(fabs(Oldres)<2.*fabs(res))
		      {
			REAL del = (Oldres*Oldres-4.*Oldres*res);
			if(del>=0.)
			  {
			    REAL l1 = (Oldres+del)/2./res;
			    REAL l2 = (Oldres-del)/2./res;
			    if(minlambda<l1 && l1<=1.)
			      lambda=l1;
			    else if(minlambda<l2 && l2<=1)
			      lambda=l2;
			    else
			      lambda = minlambda;
			  }
			else
			  lambda = Oldres/2./res;
			if(lambda<minlambda)
			  lambda=minlambda;
		      }
		    for (int var=0; var<DOF; var++)
		      {
			INDEX_PETSC(LambdaprimVecLocal,&zone,var)=
			  (1.-lambda)*INDEX_PETSC(TmpprimVecLocal,&zone,var)
			  +lambda*INDEX_PETSC(primVecLocal,&zone,var);
		      }
		  }
	      }
	    if(it==0)
	      OrigRes=DiagOldRes;
	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
				   ts->primPetscVecLambda,
				   &LambdaprimVecLocal);
	    errCode = computeResidual(ts->snes,ts->primPetscVecLambda,ts->LambdaresidualPetscVec,ts);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LambdaresidualPetscVec,
			       &LambdaresidualLocal);
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
		    REAL Lambdares = 0.;
		    REAL Oldres = 0.;
		    REAL res = 0.;
		    REAL temp = 0.;
		    for (int var=0; var<DOF; var++)
		      {
			temp = INDEX_PETSC(LambdaresidualLocal,&zone,var);
			Lambdares+=temp*temp;
			temp = INDEX_PETSC(OldresidualLocal,&zone,var);
                        Oldres+=temp*temp;
                        temp = INDEX_PETSC(residualLocal,&zone,var);
                        res+=temp*temp;
		      }
		    Lambdares=sqrt(Lambdares);
		    res=sqrt(res);
                    Oldres=sqrt(Oldres);
		    DiagLamRes += Lambdares;
		    if(res<Lambdares && res<Oldres)
		      {}//don't use line search
		    else
		      {
			if(Oldres<Lambdares && itlam>0)
			  for (int var=0; var<DOF; var++)
			    INDEX_PETSC(primVecLocal,&zone,var)=
			      INDEX_PETSC(TmpprimVecLocal,&zone,var);
			else
			  for (int var=0; var<DOF; var++)
			    INDEX_PETSC(primVecLocal,&zone,var)=
			      INDEX_PETSC(LambdaprimVecLocal,&zone,var);
		      }
		  }  
	      }
	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
				   ts->primPetscVecHalfStep,
				   &primVecLocal);
	    errCode = computeResidual(ts->snes,ts->primPetscVecHalfStep,ts->residualPetscVec,ts);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
			       &residualLocal);
	    REAL DiagFinalRes = 0.;
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
		    DiagFinalRes += res;
		  }  
	      }
	    printf("OldRes = %e; FullRes = %e; LamRes = %e; FinalRes = %e\n",DiagOldRes,DiagFullRes,DiagLamRes,DiagFinalRes);
	    if(DiagFinalRes<atol+rtol*OrigRes)
	      break;
	  }
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
    
    VecCopy(ts->primPetscVecHalfStep, ts->primPetscVec);
    for(int it=0;it<MAX_IT;it++)
      {
	VecCopy(ts->primPetscVec, ts->primPetscVecTmp);
	errCode = computeResidual(ts->snes,ts->primPetscVec,ts->OldresidualPetscVec,ts);
	errCode = SNESSolve(ts->snes, NULL, ts->primPetscVec);
	int AllPointsConverged = 0;
	REAL minlambda = 1.e-5;

	int MAX_LAMBDA_LOOPS = 3;
	DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
			   &OldresidualLocal);
	DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->primPetscVecTmp,
			   &TmpprimVecLocal);
	for(int itlam = 0; itlam<MAX_LAMBDA_LOOPS && AllPointsConverged==0;itlam++)
	  {
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
			       &residualLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LambdaresidualPetscVec,
                               &LambdaresidualLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->primPetscVec,
			   &primVecLocal);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->primPetscVecLambda,
                               &LambdaprimVecLocal);
		
	    REAL DiagOldRes = 0.;
	    REAL DiagFullRes = 0.;
	    REAL DiagLamRes = 0.;
	    AllPointsConverged = 1;
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
		    REAL Oldres = 0.;
		    REAL temp = 0.;
		    for (int var=0; var<DOF; var++)
		      {
			temp = INDEX_PETSC(OldresidualLocal,&zone,var);
			Oldres+=temp*temp;
			temp = INDEX_PETSC(residualLocal,&zone,var);
			res+=temp*temp;
		      }
		    res=sqrt(res);
		    Oldres=sqrt(Oldres);
		    DiagOldRes += Oldres;
		    DiagFullRes += res;
		    if(fabs(res)<fabs(Oldres))
		      AllPointsConverged=0;
		    REAL lambda = 1.;
		    if(fabs(Oldres)<2.*fabs(res))
		      {
			REAL del = (Oldres*Oldres-4.*Oldres*res);
			if(del>=0.)
			  {
			    REAL l1 = (Oldres+del)/2./res;
			    REAL l2 = (Oldres-del)/2./res;
			    if(minlambda<l1 && l1<=1.)
			      lambda=l1;
			    else if(minlambda<l2 && l2<=1)
			      lambda=l2;
			    else
			      lambda = minlambda;
			  }
			else
			  lambda = Oldres/2./res;
			if(lambda<minlambda)
			  lambda=minlambda;
		      }
		    for (int var=0; var<DOF; var++)
		      {
			INDEX_PETSC(LambdaprimVecLocal,&zone,var)=
			  (1.-lambda)*INDEX_PETSC(TmpprimVecLocal,&zone,var)
			  +lambda*INDEX_PETSC(primVecLocal,&zone,var);
		      }
		  }
	      }
	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
				   ts->primPetscVecLambda,
				   &LambdaprimVecLocal);
	    errCode = computeResidual(ts->snes,ts->primPetscVecLambda,ts->LambdaresidualPetscVec,ts);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LambdaresidualPetscVec,
			       &LambdaresidualLocal);
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
		    REAL Lambdares = 0.;
                    REAL Oldres = 0.;
		    REAL temp = 0.;
		    for (int var=0; var<DOF; var++)
		      {
			temp = INDEX_PETSC(LambdaresidualLocal,&zone,var);
			Lambdares+=temp*temp;
			temp = INDEX_PETSC(OldresidualLocal,&zone,var);
                        Oldres+=temp*temp;
                        temp = INDEX_PETSC(residualLocal,&zone,var);
                        res+=temp*temp;
		      }
		    Lambdares=sqrt(Lambdares);
                    res=sqrt(res);
                    Oldres=sqrt(Oldres);
                    DiagLamRes += Lambdares;
		    if(res<Lambdares && res<Oldres)
		      {}//don't use line search
		    else
                      {
			if(Oldres<Lambdares && itlam>0)
                          for (int var=0; var<DOF; var++)
                            INDEX_PETSC(primVecLocal,&zone,var)=
                              INDEX_PETSC(TmpprimVecLocal,&zone,var);
			else
			  for (int var=0; var<DOF; var++)
			    INDEX_PETSC(primVecLocal,&zone,var)=
			      INDEX_PETSC(LambdaprimVecLocal,&zone,var);
                      }
		  }  
	      }
	    if(it==0)
	      OrigRes=DiagOldRes;
	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
				   ts->primPetscVec,
				   &primVecLocal);
	    errCode = computeResidual(ts->snes,ts->primPetscVec,ts->residualPetscVec,ts);
	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
			       &residualLocal);
	    REAL DiagFinalRes = 0.;
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
		    DiagFinalRes += res;
		  }  
	      }
	    printf("OldRes = %e; FullRes = %e; LamRes = %e; FinalRes = %e\n",DiagOldRes,DiagFullRes,DiagLamRes,DiagFinalRes);
	    if(DiagFinalRes<atol+rtol*OrigRes)
	      break;
	  }
      }
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
  VecDestroy(&ts->primPetscVecTmp);
  VecDestroy(&ts->primPetscVecLambda);
  VecDestroy(&ts->divFluxPetscVecOld);
  VecDestroy(&ts->conservedVarsPetscVecOld);
  VecDestroy(&ts->sourceTermsPetscVecOld);
  VecDestroy(&ts->residualPetscVec);
  VecDestroy(&ts->OldresidualPetscVec);
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
