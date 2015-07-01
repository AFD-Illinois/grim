//#include "timestepper.h"
//
//
///* Explicit time stepping:
// * -----------------------
// * The time stepping is done in two stages:
// *
// * 1) Get prim at t=n+1/2 using:
// *    
// *    (U^(n+1/2) - U^n)/(.5*dt) + grad(F)^(n-1) + sources^(n-1) = 0
// *    
// * 2) Now using prim at t=n+1/2, compute fluxes and sources and go from t=n to
// *    t=n+1:
// *
// *    (U^(n+1) - U^n)/dt + grad(F)^(n+1/2) + sources^(n+1/2) = 0
// *
// * IMEX time stepping:
// * -------------------
// * The time stepping is done in two stages:
// *
// * 1) Get prim at t=n+1/2 using:
// *    
// *    (U^(n+1/2) - U^n)/(.5*dt) + grad(F)^(n-1) + sources^(n-1) = 0
// *    
// * 2) Now using prim at t=n+1/2, compute fluxes but compute sources as a
// *    combination of sources at t=n and t=n+1 and go from t=n to t=n+1:
// *
// *    (U^(n+1) - U^n)/dt + grad(F)^(n+1/2) + 0.5*(sources^n + sources^(n+1)) = 0
// * 
// * Implicit time stepping:
// * -----------------------
// * Go directly from t=n to t=n+1 using
// * 
// *    (U^(n+1) - U^n)/dt + 0.5*(grad(F)^n + grad(F)^(n+1))
// *                       + 0.5*(sources^n + sources^(n+1)) = 0
// */
//void timeStep(struct timeStepper ts[ARRAY_ARGS 1])
//{
//  printf("Start time step \n");
//  PetscErrorCode errCode;
//
//  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
//
//    /* First go from t=n to t=n+1/2 */
//    ts->computeOldSourceTermsAndOldDivOfFluxes = 1;
//    ts->computeDivOfFluxAtTimeN = 1;
//    ts->computeDivOfFluxAtTimeNPlusHalf = 0;
//    ts->computeSourceTermsAtTimeN = 1;
//    ts->computeSourceTermsAtTimeNPlusHalf = 0;
//
//    REAL dt = ts->dt;
//    ts->dt = dt/2.;
//
//    int MAX_IT = 3;
//    int MAX_LAMBDA_LOOPS = 10;
//    double ALPHA_LS = 1.e-4;
//    double LAMBDA_MIN = 1.e-12;
//    REAL atol = 1.e-7;
//    REAL rtol = 1.e-5;
//
//    ARRAY(OldresidualLocal);
//    ARRAY(LambdaresidualLocal);
//    ARRAY(LastStepresidualLocal);
//    ARRAY(residualLocal);
//    ARRAY(HalfStepprimVecLocal);
//    ARRAY(OldprimVecLocal);
//    ARRAY(LastStepprimVecLocal);
//    ARRAY(primVecLocal);
//    ARRAY(LambdaprimVecLocal);
//    
//    
//    int Npoints = ts->X1Size*ts->X2Size;
//    int UseCustomLineSearch = 1;
//    REAL mLambda[Npoints];
//    int mConverged[Npoints];
//
//    if(UseCustomLineSearch==1)
//      {
//	errCode = computeResidual(ts->snes,ts->primPetscVecOld,ts->OldresidualPetscVec,ts);
//	VecCopy(ts->primPetscVecOld, ts->primPetscVecHalfStep);
//	VecCopy(ts->OldresidualPetscVec, ts->residualPetscVec);
//	DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
//			   &OldresidualLocal);
//	DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecOld,
//			   &OldprimVecLocal);
//	//Newton linesearch method
//	for(int it=0;it<MAX_IT;it++)
//	  {
//	    VecCopy(ts->residualPetscVec,ts->LastStepresidualPetscVec);
//	    VecCopy(ts->primPetscVecHalfStep, ts->primPetscVecLastStep);
//	    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecLastStep,
//			       &LastStepprimVecLocal);
//	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LastStepresidualPetscVec,
//			       &LastStepresidualLocal);
//	    //1) Standard petsc inversion. Should use basic linesearch, and snes_max_it=1
//	    errCode = SNESSolve(ts->snes, NULL, ts->primPetscVecHalfStep);
//	    int AllPointsConverged = 0;
//	    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
//			       &primVecLocal);
//	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
//			       &residualLocal);
//	    VecCopy(ts->residualPetscVec,ts->LambdaresidualPetscVec);
//	    VecCopy(ts->primPetscVecHalfStep, ts->primPetscVecLambda);
//	    //Linesearch within each step
//	    for(int itlam = 0; itlam<MAX_LAMBDA_LOOPS && AllPointsConverged==0;itlam++)
//	      {
//		DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecLambda,
//				   &LambdaprimVecLocal);
//		DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LambdaresidualPetscVec,
//				   &LambdaresidualLocal);
//		//Compute L2 norm of residual pointwise
//		AllPointsConverged = 1;
//		int NptsNotConv = 0;
//		REAL DiagOldRes = 0.;
//		REAL DiagNewRes = 0.;
//                #if (USE_OPENMP)
//                  #pragma omp parallel for
//                #endif
//		LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
//		  {
//		    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
//		      {
//			int idx = iTile*ts->X2Size*TILE_SIZE_X1
//			  +jTile*TILE_SIZE_X1*TILE_SIZE_X2+
//			  iInTile*TILE_SIZE_X2+jInTile;
//			if(itlam==0)
//			  {
//			    mLambda[idx]=1.;
//			    mConverged[idx]=0;
//			  }
//			struct gridZone zone;
//			setGridZone(iTile, jTile,
//				    iInTile, jInTile,
//				    ts->X1Start, ts->X2Start,
//				    ts->X1Size, ts->X2Size,
//				    &zone);
//			REAL Oldres = 0.;
//			REAL Lambdares = 0.;
//			REAL Firstres = 0.;
//			REAL temp = 0.;
//			for (int var=0; var<DOF; var++)
//			  {
//			    temp = INDEX_PETSC(LastStepresidualLocal,&zone,var);
//			    Oldres+=temp*temp;
//			    temp = INDEX_PETSC(LambdaresidualLocal,&zone,var);
//			    Lambdares+=temp*temp; 
//			    temp = INDEX_PETSC(OldresidualLocal,&zone,var);
//			    Firstres+=temp*temp;
//			  }
//			DiagOldRes += Oldres;
//			DiagNewRes += Lambdares;
//			//Stop if we already converged (we went this far just
//			//to include all points in the residual)
//			if(mConverged[idx]==1)
//			  continue;
//			if(Lambdares<Oldres*(1.-mLambda[idx]*ALPHA_LS))
//			  {
//			    //Acceptable end to the linesearch
//			    mConverged[idx]=1;
//			  }
//			else
//			  {
//			    //if the residual was already small, just accept it...
//			    if(sqrt(Oldres)<atol+rtol*sqrt(Firstres))
//			      {
//				for (int var=0; var<DOF; var++)
//				  INDEX_PETSC(LambdaprimVecLocal,&zone,var)= 
//				    INDEX_PETSC(LastStepprimVecLocal,&zone,var);
//				mConverged[idx]=1;
//				mLambda[idx]=0.;
//				DiagNewRes+=Oldres-Lambdares;
//			      }
//			    else
//			      {
//				//backtrack
//				REAL lam = mLambda[idx];
//				REAL lamupd = Oldres*lam*lam/(Lambdares+(2.*lam-1.)*Oldres);
//				if(lamupd>0.5*lam)
//				  lamupd=0.5*lam;
//				if(lamupd<0.1*lam)
//				  lamupd=0.1*lam;
//				mLambda[idx]=lamupd;
//				if(mLambda[idx]<LAMBDA_MIN)
//				  mLambda[idx]=LAMBDA_MIN;
//				for (int var=0; var<DOF; var++)
//				  INDEX_PETSC(LambdaprimVecLocal,&zone,var)=
//				    (1.-mLambda[idx])*INDEX_PETSC(LastStepprimVecLocal,&zone,var)+
//				    mLambda[idx]*INDEX_PETSC(primVecLocal,&zone,var);
//				AllPointsConverged=0;
//				NptsNotConv++;
//			      }
//			  }
//		      }
//		  }
//		DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
//				       ts->primPetscVecLambda,
//				       &LambdaprimVecLocal);
//		DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, 
//				       ts->LambdaresidualPetscVec,
//				       &LambdaresidualLocal);
//		DiagOldRes = sqrt(DiagOldRes);
//		DiagNewRes = sqrt(DiagNewRes);
//		if(itlam==0)
//		  printf("Step %i - Initial residual = %e \n",it,DiagOldRes);
//		printf("New residual = %e. ",DiagNewRes);
//		printf("%i of %i points have not converged.\n",NptsNotConv,Npoints);
//		errCode = computeResidual(ts->snes,ts->primPetscVecLambda,ts->LambdaresidualPetscVec,ts);
//	      }
//	    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
//				   ts->primPetscVecHalfStep,
//				   &primVecLocal);
//	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, 
//				   ts->residualPetscVec,
//				   &residualLocal);
//	    VecCopy(ts->LambdaresidualPetscVec,ts->residualPetscVec);
//	    VecCopy(ts->primPetscVecLambda, ts->primPetscVecHalfStep);
//	    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
//				   ts->primPetscVecLastStep,
//				   &LastStepprimVecLocal);
//	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, 
//				   ts->LastStepresidualPetscVec,
//				   &LastStepresidualLocal);
//	  }
//	DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
//			       &OldresidualLocal);
//	DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecOld,
//			       &OldprimVecLocal);
//      }
//    else
//      {
//	VecCopy(ts->primPetscVecOld, ts->primPetscVecHalfStep);
//	errCode = SNESSolve(ts->snes, NULL, ts->primPetscVecHalfStep);
//      }
//
//    /* Problem dependent half step diagnostics */
//    halfStepDiagnostics(ts);
//
//    /* Current state:
//    * ts->primPetscVecHalfStep has the primitive variables at t = n+1/2 
//    * ts->primPetscVecOld has the primitive variables at t = n
//    */
//
//    /* Now go from t=n to t=n+1 using
//    * 1) EXPLICIT time stepping:
//    * fluxes and sources computed using primitive variables at t=n+1/2.
//    * 
//    * 2) IMEX time stepping
//    * fluxes computed using primitive variables at t=n+1/2 whereas sources
//    * computed using 0.5*(Sources^n + Sources^(n-1))*/
//    ts->computeOldSourceTermsAndOldDivOfFluxes = 1;
//    ts->computeDivOfFluxAtTimeN = 0;
//    ts->computeDivOfFluxAtTimeNPlusHalf = 1;
//
//    #if (TIME_STEPPING==EXPLICIT)
//      ts->computeSourceTermsAtTimeN = 0;
//      ts->computeSourceTermsAtTimeNPlusHalf = 1;
//    #elif (TIME_STEPPING==IMEX)
//      ts->computeSourceTermsAtTimeN = 0;
//      ts->computeSourceTermsAtTimeNPlusHalf = 0;
//    #endif
//
//    ts->dt = dt;
//
//    if(UseCustomLineSearch==1)
//      {
//	errCode = computeResidual(ts->snes,ts->primPetscVecHalfStep,ts->OldresidualPetscVec,ts);
//	VecCopy(ts->primPetscVecHalfStep, ts->primPetscVec);
//	VecCopy(ts->OldresidualPetscVec, ts->residualPetscVec);
//	DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
//			   &OldresidualLocal);
//	DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
//			   &HalfStepprimVecLocal);
//	for(int it=0;it<MAX_IT;it++)
//	  {
//	    VecCopy(ts->residualPetscVec,ts->LastStepresidualPetscVec);
//	    VecCopy(ts->primPetscVec, ts->primPetscVecLastStep);
//	    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecLastStep,
//			       &LastStepprimVecLocal);
//	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LastStepresidualPetscVec,
//			       &LastStepresidualLocal);
//	    errCode = SNESSolve(ts->snes, NULL, ts->primPetscVec);
//	    int AllPointsConverged = 0;
//	    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVec,
//			       &primVecLocal);
//	    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
//			       &residualLocal);
//	    VecCopy(ts->residualPetscVec,ts->LambdaresidualPetscVec);
//	    VecCopy(ts->primPetscVec, ts->primPetscVecLambda);
//	    //Linesearch within each step
//	    for(int itlam = 0; itlam<MAX_LAMBDA_LOOPS && AllPointsConverged==0;itlam++)
//	      {
//		DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecLambda,
//				   &LambdaprimVecLocal);
//		DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->LambdaresidualPetscVec,
//				   &LambdaresidualLocal);
//		//Compute L2 norm of residual pointwise                                                                                                                         
//		AllPointsConverged = 1;
//		int NptsNotConv = 0;
//		REAL DiagOldRes = 0.;
//		REAL DiagNewRes = 0.;
//                #if (USE_OPENMP)
//                  #pragma omp parallel for
//                #endif
//		LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
//		  {
//		    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
//		      {
//			int idx = iTile*ts->X2Size*TILE_SIZE_X1
//			  +jTile*TILE_SIZE_X1*TILE_SIZE_X2+
//			  iInTile*TILE_SIZE_X2+jInTile;
//			if(itlam==0)
//			  {
//			    mLambda[idx]=1.;
//			    mConverged[idx]=0;
//			  }
//			struct gridZone zone;
//			setGridZone(iTile, jTile,
//				    iInTile, jInTile,
//				    ts->X1Start, ts->X2Start,
//				    ts->X1Size, ts->X2Size,
//				    &zone);
//			REAL Oldres = 0.;
//			REAL Lambdares = 0.;
//			REAL Firstres = 0.;
//			REAL temp = 0.;
//			for (int var=0; var<DOF; var++)
//			  {
//			    temp = INDEX_PETSC(OldresidualLocal,&zone,var);
//			    Firstres+=temp*temp;
//			    temp = INDEX_PETSC(LastStepresidualLocal,&zone,var);
//			    Oldres+=temp*temp;
//			    temp = INDEX_PETSC(LambdaresidualLocal,&zone,var);
//			    Lambdares+=temp*temp; 
//			  }
//			DiagOldRes += Oldres;
//			DiagNewRes += Lambdares;
//			//Stop if we already converged (we went this far just
//			//to include all points in the residual)
//			if(mConverged[idx]==1)
//			  continue;
//			if(Lambdares<Oldres*(1.-mLambda[idx]*ALPHA_LS))
//			  {
//			    //Acceptable end to the linesearch
//			    mConverged[idx]=1;
//			  }
//			else
//			  {
//			    //if the residual was already small, just accept it...
//			    if(sqrt(Oldres)<atol+rtol*sqrt(Firstres))
//			      {
//				for (int var=0; var<DOF; var++)
//				  INDEX_PETSC(LambdaprimVecLocal,&zone,var)= 
//				    INDEX_PETSC(LastStepprimVecLocal,&zone,var);
//				mConverged[idx]=1;
//				mLambda[idx]=0.;
//				DiagNewRes+=Oldres-Lambdares;
//			      }
//			    else
//			      {
//				//backtrack
//				REAL lam = mLambda[idx];
//				REAL lamupd = Oldres*lam*lam/(Lambdares+(2.*lam-1.)*Oldres);
//				if(lamupd>0.5*lam)
//                                  lamupd=0.5*lam;
//                                if(lamupd<0.1*lam)
//                                  lamupd=0.1*lam;
//                                mLambda[idx]=lamupd;
//				if(mLambda[idx]<LAMBDA_MIN)
//				  mLambda[idx]=LAMBDA_MIN;
//				for (int var=0; var<DOF; var++)
//				  INDEX_PETSC(LambdaprimVecLocal,&zone,var)=
//				    (1.-mLambda[idx])*INDEX_PETSC(LastStepprimVecLocal,&zone,var)+
//				    mLambda[idx]*INDEX_PETSC(primVecLocal,&zone,var);
//				AllPointsConverged=0;
//				NptsNotConv++;
//			      }
//			  }
//		      }
//		  }
//		DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
//				       ts->primPetscVecLambda,
//				       &LambdaprimVecLocal);
//		DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
//				       ts->LambdaresidualPetscVec,
//				       &LambdaresidualLocal);
//		DiagOldRes = sqrt(DiagOldRes);
//		DiagNewRes = sqrt(DiagNewRes);
//		if(itlam==0 && it==0)
//		  printf("Initial residual = %e\n",DiagOldRes);
//		if(itlam==0)
//		  printf("Step %i \n",it);
//		printf("New residual = %e. ",DiagNewRes);
//		printf("%i of %i points have not converged.\n",NptsNotConv,Npoints);
//		errCode = computeResidual(ts->snes,ts->primPetscVecLambda,ts->LambdaresidualPetscVec,ts);
//	      }
//	    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
//				   ts->primPetscVec,
//				   &primVecLocal);
//	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, 
//				   ts->residualPetscVec,
//				   &residualLocal);
//	    VecCopy(ts->LambdaresidualPetscVec,ts->residualPetscVec);
//	    VecCopy(ts->primPetscVecLambda, ts->primPetscVec);
//	    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
//				   ts->primPetscVecLastStep,
//				   &LastStepprimVecLocal);
//	    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
//				   ts->LastStepresidualPetscVec,
//				   &LastStepresidualLocal);
//	  }
//	DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, ts->OldresidualPetscVec,
//			       &OldresidualLocal);
//	DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
//			       &HalfStepprimVecLocal);
//      }    
//    else
//      {
//	VecCopy(ts->primPetscVecHalfStep, ts->primPetscVec);
//	errCode = SNESSolve(ts->snes, NULL, ts->primPetscVec);
//      }
//    //Diagnostics
//    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVec,
//		       &primVecLocal);
//    DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, ts->residualPetscVec,
//		       &residualLocal);   
//    int NptsResMag[16];
//    for(int i=0;i<16;i++)
//      NptsResMag[i]=0;
//    #if (USE_OPENMP)
//      #pragma omp parallel for
//    #endif
//    LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
//      {
//	LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
//	  {
//	    struct gridZone zone;
//	    setGridZone(iTile, jTile,
//			iInTile, jInTile,
//			ts->X1Start, ts->X2Start,
//			ts->X1Size, ts->X2Size,
//			&zone);
//	    REAL res = 0.;
//	    REAL temp = 0.;
//	    for (int var=0; var<DOF; var++)
//	      {
//		temp = INDEX_PETSC(residualLocal,&zone,var);
//		res+=temp*temp;
//	      }
//	    res=sqrt(res);
//	    res=log10(res+1.e-15);
//	    for(int i=0;i<16;i++)
//	      if(res>-i)
//		{
//		  NptsResMag[i]++;
//		  break;
//		}
//	  }
//      }
//    printf("Distribution of residuals : ");
//    for(int i=0;i<16;i++)
//      printf("%i ;", NptsResMag[i]);
//    printf("\n");
//    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
//			   ts->primPetscVec,
//			   &primVecLocal);
//    DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
//			   ts->residualPetscVec,
//			   &residualLocal);    
//
//  #elif (TIME_STEPPING==IMPLICIT)
//
//    ts->computeOldSourceTermsAndOldDivOfFluxes = 1;
//    ts->computeDivOfFluxAtTimeN = 1;
//    ts->computeDivOfFluxAtTimeNPlusHalf = 0;
//    ts->computeSourceTermsAtTimeN = 1;
//    ts->computeSourceTermsAtTimeNPlusHalf = 0;
//
//    VecCopy(ts->primPetscVecOld, ts->primPetscVec);
//    errCode = SNESSolve(ts->snes, NULL, ts->primPetscVec);
//
//  #endif
//
//  ts->t = ts->t + ts->dt;
//  ts->timeStepCounter++;
//  PetscPrintf(PETSC_COMM_WORLD,
//              "\nCompleted step %d, current time = %.5f, dt = %.5f\n\n",
//              ts->timeStepCounter, ts->t, ts->dt);
//
//  /* Problem dependent full step diagnostics */
//  fullStepDiagnostics(ts);
//
//  /* Copy solved variables to old variables */
//  VecCopy(ts->primPetscVec, ts->primPetscVecOld);
//
//  diagnostics(ts);
//
//  REAL newDt;
//  #if (COMPUTE_DIM==1)
//
//    VecStrideMin(ts->dtPetscVec, 0, NULL, &newDt);
//
//  #elif (COMPUTE_DIM==2)
//    REAL dtX1, dtX2;
//    VecStrideMin(ts->dtPetscVec, 0, NULL, &dtX1);
//    VecStrideMin(ts->dtPetscVec, 1, NULL, &dtX2);
//
//    newDt = 1./( (1./dtX1) + (1./dtX2) );
//  #endif
//
//  if (newDt > MAX_DT_INCREMENT*ts->dt)
//  {
//    newDt = MAX_DT_INCREMENT*ts->dt;
//  }
//
//  ts->dt = newDt;
//}
