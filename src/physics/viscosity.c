#include "physics.h"

void addViscositySourceTermsToResidual
(
  const REAL primTile[ARRAY_ARGS TILE_SIZE],
  ARRAY(primGlobal), ARRAY(primHalfStepGlobal), ARRAY(primOldGlobal),
  ARRAY(connectionGlobal),
  ARRAY(graduConVisGlobal), 
  ARRAY(graduConHigherOrderTerm1VisGlobal),
  REAL dt,
  int computeOldSourceTermsAndOldDivOfFluxes,
  int computeDivOfFluxAtTimeN,
  int computeDivOfFluxAtTimeNPlusHalf,
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  ARRAY(residualGlobal)
)
{
#if (VISCOSITY)
  LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
  {
    struct gridZone zoneCenter;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zoneCenter);
    
    /* Values needed to compute spatial derivative.  Depending on where the
     * spatial derivatives are to be computed (t=n or t=n+1/2), the values of
     * primTile is set (already done, outside this routine) */

    if (computeOldSourceTermsAndOldDivOfFluxes)
    {
      REAL graduConVis[COMPUTE_DIM*NDIM];
      REAL graduConHigherOrderTerm1Vis[COMPUTE_DIM];
      computeViscositySpatialGradientTerms
        (primTile, 
         iTile, jTile, 
         iInTile, jInTile,
         X1Start, X2Start,
         X1Size, X2Size, 
         graduConVis, 
         graduConHigherOrderTerm1Vis
        );
      
      for (int mu=0; mu<NDIM; mu++)
      {
        INDEX_PETSC(graduConVisGlobal, &zoneCenter, mu) = graduConVis[mu];
      }

      INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 0) = 
        graduConHigherOrderTerm1Vis[0];

      #if (COMPUTE_DIM==2)
        for (int mu=0; mu<NDIM; mu++)
        {
          INDEX_PETSC(graduConVisGlobal, &zoneCenter, mu+NDIM) = graduConVis[mu+NDIM];
        }

        INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 1) = 
          graduConHigherOrderTerm1Vis[1];

      #endif
    } 
    else
    {
      #if (TIME_STEPPING==IMPLICIT)
        REAL graduConVis[COMPUTE_DIM*NDIM];
        REAL graduConHigherOrderTerm1Vis[COMPUTE_DIM];
        computeViscositySpatialGradientTerms
          (primTile, 
           iTile, jTile, 
           iInTile, jInTile,
           X1Start, X2Start,
           X1Size, X2Size, 
           graduConVis,
           graduConHigherOrderTerm1Vis,
          );
      #endif

      /* Setup all the data structures needed 
       * Note: EXPLICIT and IMEX time stepping algorithms have a half step
       * involved whereas the IMPLICIT time step is a single step algorithm 
       * 
       * Legend: elem       --> fluid element at the (n+1)th time step 
       *
       *         elemCenter --> fluid element at the nth time step when in the
       *                        first  half of EXPLICIT/IMEX time step and
       *                        (n+1/2) time step when in the second half the
       *                        the EXPLICIT/IMEX time step. SHOULD NOT BE USED
       *                        WITH IMPLICIT time stepping since it is not set.
       *
       *         elemOld    --> fluid element at the nth time step.
       *
       *         zoneCenter --> zone at (iInTile, jInTile).
       *
       *         geomCenter --> geometry at the (iInTile, jInTile) zone.
       */
      struct fluidElement elem, elemCenter, elemOld;
      struct geometry geomCenter;

      REAL XCoords[NDIM];
      getXCoords(&zoneCenter, CENTER, XCoords);
      setGeometry(XCoords, &geomCenter);

      setFluidElement(&INDEX_PETSC(primGlobal, &zoneCenter, 0), &geomCenter,
                      &elem);

      setFluidElement(&INDEX_PETSC(primOldGlobal, &zoneCenter, 0), &geomCenter,
                      &elemOld);

      if (computeDivOfFluxAtTimeN)
      {
        setFluidElement(&INDEX_PETSC(primOldGlobal, &zoneCenter, 0), 
                        &geomCenter, &elemCenter);
      } 
      else if (computeDivOfFluxAtTimeNPlusHalf)
      {
        setFluidElement(&INDEX_PETSC(primHalfStepGlobal, &zoneCenter, 0),
                        &geomCenter, &elemCenter);
      }

      //Temperature information
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)
        REAL TCenter =   (ADIABATIC_INDEX-1.)
                       * elemCenter.primVars[UU]
                       / elemCenter.primVars[RHO];
      #elif (TIME_STEPPING == IMPLICIT)
	REAL TCenter = (ADIABATIC_INDEX-1.)
	  * elem.primVars[UU]
	  / elem.primVars[RHO];
      #endif
      if(TCenter<1.e-12)
	TCenter = 1.e-12;

	//Compute beta at Center
	//where beta = tauVis/2/eta
      #if (HIGHORDERTERMS_VISCOSITY)
        #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)
          REAL betaCenter = elemCenter.tauVis/(elemCenter.eta*2.);
        #elif (TIME_STEPPING == IMPLICIT)
          REAL betaCenter = 0.5*(  (elem.tauVis/(elem.eta*2.) )
				 + (elemOld.tauVis/(elemOld.eta*2.) )
				 );
        #endif
      #endif

      REAL divUCoeff = 1.;
      #if (HIGHORDERTERMS_VISCOSITY)
        divUCoeff = .5;
      #endif

      /* Higher order term 1 
       * Psi/2*(sqrt(g)*u^\mu)_{,\mu} */
      REAL dHigherOrderTerm1[NDIM];
      
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)

        dHigherOrderTerm1[0] = 
          divUCoeff*elemCenter.primVars[PSI] 
        * sqrt(-geomCenter.gDet) 
        * (elem.uCon[0] - elemOld.uCon[0])/dt;
      

        dHigherOrderTerm1[1] = 
          divUCoeff*elemCenter.primVars[PSI]
        * INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 0);

        dHigherOrderTerm1[2] = 0.; dHigherOrderTerm1[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm1[2] = 
            divUCoeff*elemCenter.primVars[PSI]
          * INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 1);
        #endif

      #elif (TIME_STEPPING == IMPLICIT)

        dHigherOrderTerm1[0] = 
          0.5 * divUCoeff * (elem.primVars[PSI] + elemOld.primVars[PSI])
              * sqrt(-geomCenter.gDet) 
              * (elem.uCon[0] - elemOld.uCon[0])/dt;
      
        dHigherOrderTerm1[1] = 
          0.5* divUCoeff * (  (  elemOld.primVars[PSI]
                  * INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 0)
                 )
               + 
                 (  elem.primVars[PSI]
                  * graduConHigherOrderTerm1Vis[0]
                 )
              );

        dHigherOrderTerm1[2] = 0.; dHigherOrderTerm1[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm1[2] = 
	    0.5* divUCoeff * (  (  elemOld.primVars[PSI]
                  * INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 1)
                 )
               + 
                (  elem.primVars[PSI]
                  * graduConHigherOrderTerm1Vis[1]
                 )
	      );
        #endif

      #endif /* TIME_STEPPING options for higherOrderTerm1 */

      REAL higherOrderTerm1 = 0.;
      for (int mu=0; mu<NDIM; mu++)
        higherOrderTerm1 += dHigherOrderTerm1[mu];
      
      // We now compute the target pressure anisotropy.
      // dP_0 = 3*eta*b^mu*b^nu u_{mu;nu} - eta * u^{mu}_{;mu}
      REAL TargetPsi = 0.;
      REAL g = sqrt(-geomCenter.gDet);
      REAL norm = g;

      //Coordinate derivatives
      {
      #if (TIME_STEPPING == IMEX || TIME_STEPPING == EXPLICIT)
	REAL bSqr, bCovCenter[NDIM];
        conToCov(elemCenter.bCon, &geomCenter, bCovCenter);
        bSqr = getbSqr(&elemCenter, &geomCenter);
        for (int alpha=0; alpha<NDIM; alpha++)
	  {
	    TargetPsi+=bCovCenter[alpha]*elemCenter.bCon[0]/bSqr*3.*elem.eta
	      * (elem.uCon[alpha] - elemOld.uCon[alpha])/dt;
	    TargetPsi+=bCovCenter[alpha]*elemCenter.bCon[1]/bSqr*3.*elem.eta
	      * INDEX_PETSC(graduConVisGlobal, &zoneCenter, alpha);
            #if (COMPUTE_DIM==2)
	      TargetPsi+=bCovCenter[alpha]*elemCenter.bCon[2]/bSqr*3.*elem.eta
		* INDEX_PETSC(graduConVisGlobal, &zoneCenter, alpha+NDIM);
            #endif
	  }
	TargetPsi -= elem.eta
	  * (elem.uCon[0] - elemOld.uCon[0])/dt;
	TargetPsi -= elem.eta/sqrt(-geomCenter.gDet)
	  *INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 0);
        #if (COMPUTE_DIM==2)
          TargetPsi -= elem.eta/sqrt(-geomCenter.gDet)
	    *INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 1);
        #endif
      #elif (TIME_STEPPING == IMPLICIT)
	REAL bSqr, bCov[NDIM];
	REAL bSqrOld, bCovOld[NDIM];
	conToCov(elem.bCon, &geomCenter, bCov);
	conToCov(elemOld.bCon, &geomCenter, bCovOld);
	bSqr = getbSqr(&elem, &geomCenter);
	bSqrOld = getbSqr(&elemOld, &geomCenter);
	for (int alpha=0; alpha<NDIM; alpha++)
	  {
	    TargetPsi+=1.5*(bCov[alpha]*elem.bCon[0]*elem.eta/bSqr
			    +bCovOld[alpha]*elemOld.bCon[0]*elemOld.eta/bSqrOld)
	      * (elem.uCon[alpha] - elemOld.uCon[alpha])/dt;
	    TargetPsi+=1.5*(bCov[alpha]*elem.bCon[1]*elem.eta*graduConVis[alpha]/bSqr
			    +bCovOld[alpha]*elemOld.bCon[1]*elemOld.eta/bSqrOld
			    *INDEX_PETSC(graduConVisGlobal, &zoneCenter, alpha));
            #if (COMPUTE_DIM==2)
              TargetPsi+=1.5*(bCov[alpha]*elem.bCon[2]*elem.eta*graduConVis[alpha]/bSqr
			      +bCovOld[alpha]*elemOld.bCon[2]*elemOld.eta/bSqrOld
			      *INDEX_PETSC(graduConVisGlobal, &zoneCenter, alpha));
            #endif
	  }
	TargetPsi -= elem.eta
          * (elem.uCon[0] - elemOld.uCon[0])/dt;
	TargetPsi -= elem.eta/sqrt(-geomCenter.gDet)
          *0.5*(INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 0)+
		graduConHigherOrderTerm1Vis[0]);
        #if (COMPUTE_DIM==2)
	  TargetPsi -= elem.eta/sqrt(-geomCenter.gDet)
	    *0.5*(INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 1)+
		  graduConHigherOrderTerm1Vis[1]);
        #endif
      #endif
      }
      //Add connections
      {
      #if (TIME_STEPPING == IMEX || TIME_STEPPING == IMPLICIT)
      	REAL bSqr, bCov[NDIM];
	REAL bSqrOld, bCovOld[NDIM];
	conToCov(elem.bCon, &geomCenter, bCov);
	conToCov(elemOld.bCon, &geomCenter, bCovOld);
	bSqr = getbSqr(&elem, &geomCenter);
	bSqrOld = getbSqr(&elemOld, &geomCenter);
	REAL eOverB2 = elem.eta/bSqr;
	REAL OldeOverB2 = elemOld.eta/bSqrOld;
	for(int mu=0; mu<NDIM; mu++)
	  {
	    REAL bn = bCov[mu]*eOverB2;
	    REAL Oldbn = bCovOld[mu]*OldeOverB2;
	    for (int alpha=0; alpha<NDIM; alpha++)
	      {
		for (int gamma=0; gamma<NDIM; gamma++)
		  {
		    TargetPsi+=
		      1.5*INDEX_PETSC(connectionGlobal, &zoneCenter, GAMMA_UP_DOWN_DOWN(mu, alpha, gamma))
		      *(
			elem.bCon[alpha]*elem.uCon[gamma]*bn
			+
			elemOld.bCon[alpha]*elemOld.uCon[gamma]*Oldbn
			);
		  }
	      }
	  }
      #elif (TIME_STEPPING == EXPLICIT)
	REAL bSqr, bCov[NDIM];
        conToCov(elemCenter.bCon, &geomCenter, bCov);
        bSqr = getbSqr(&elemCenter, &geomCenter);
        for (int alpha=0; alpha<NDIM; alpha++)
          {
            for (int gamma=0; gamma<NDIM; gamma++)
	      {
		for(int mu=0; mu<NDIM; mu++)
                  {
                    TargetPsi+=
                      3.*INDEX_PETSC(connectionGlobal, &zoneCenter, GAMMA_UP_DOWN_DOWN(mu, alpha, gamma))
                      *(elemCenter.bCon[alpha]*bCov[mu]*elemCenter.uCon[gamma]/bSqr*elemCenter.eta);
                  }
	      }
	  }
      #endif	
      }
      //Now transform from dP to psi = dP * sqrt(beta/T)

      #if (HIGHORDERTERMS_VISCOSITY)
        TargetPsi*=sqrt(betaCenter/TCenter);
      #endif

      //Put the residual together
      #if (TIME_STEPPING == EXPLICIT)

      INDEX_PETSC(residualGlobal, &zoneCenter, PSI)*=elem.tauVis;

        INDEX_PETSC(residualGlobal, &zoneCenter, PSI) += 
          (- higherOrderTerm1*elem.tauVis
           + g*( elemCenter.primVars[PSI]
                 - TargetPsi
               )
          )/norm;

      #elif (TIME_STEPPING == IMEX || TIME_STEPPING == IMPLICIT)
  
	REAL flux = INDEX_PETSC(residualGlobal, &zoneCenter, PSI);
	INDEX_PETSC(residualGlobal, &zoneCenter, PSI) += 
	  (- higherOrderTerm1 //*elem.tauVis
	   + g/elemCenter.tauVis
	   *(0.5*(elem.primVars[PSI]+elemOld.primVars[PSI]) 
	   //+ g*( 0.5*(elem.primVars[PSI] + elemOld.primVars[PSI])
		 - TargetPsi
		 )
	   )/norm;
	INDEX_PETSC(residualGlobal, &zoneCenter, PSI)*=elemCenter.tauVis;

	if(0)
	  {
	    if(iTile == 2 && jTile == 4 && iInTile == 3 && jInTile == 5)
	      {
		REAL xCoords[NDIM];
		XTox(geomCenter.XCoords, xCoords);
		REAL bCov[NDIM], bSqr, uCov[NDIM];
		bSqr = getbSqr(&elem, &geomCenter);
		conToCov(elem.uCon, &geomCenter, uCov);
		conToCov(elem.bCon, &geomCenter, bCov);
		
		printf("Vars = %e; %e; %e,%e,%e; %e,%e,%e; %e\n",
		       elem.primVars[RHO],
		       elem.primVars[UU],
		       elem.primVars[U1],
		       elem.primVars[U2],
		       elem.primVars[U3],
		       elem.primVars[B1],
		       elem.primVars[B2],
		       elem.primVars[B3],
		       elem.primVars[PSI]);
		printf("Gamma = %e; uCon[0] = %e; uCov[1]=%e; bSqr = %e; Tau = %e\n",
		       elem.gamma,elem.uCon[0],uCov[1],bSqr,elem.tauVis);
		printf("Residual = %e; Flux = %e; Target = %e; HO = %e\n",
		       INDEX_PETSC(residualGlobal, &zoneCenter, PSI),
		       flux,
		       TargetPsi,
		       - higherOrderTerm1 *elem.tauVis/norm);
	      }
	  }

      #endif

      //Renormalize residual, so that it scales like dP
      #if (HIGHORDERTERMS_VISCOSITY)
	INDEX_PETSC(residualGlobal, &zoneCenter, PSI)*=sqrt(TCenter/betaCenter);
      #endif

    }
  }
#endif /* VISCOSITY */
}

void computeViscositySpatialGradientTerms
(
  const REAL primTile[ARRAY_ARGS TILE_SIZE],
  const int iTile, const int jTile,
  const int iInTile, const int jInTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  REAL graduConVis[COMPUTE_DIM*NDIM],
  REAL graduConHigherOrderTerm1Vis[COMPUTE_DIM]
)
{
#if (VISCOSITY)
  REAL XCoords[NDIM];
  struct gridZone zone, zoneCenter;
  struct geometry geom;
  struct fluidElement elem;
  REAL uConLeft[NDIM], uConCenter[NDIM], uConRight[NDIM];
  REAL gLeft, gCenter, gRight;

  setGridZone(iTile, jTile,
              iInTile, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zoneCenter);

  /* uConCenter */
  getXCoords(&zoneCenter, CENTER, XCoords);
  setGeometry(XCoords, &geom); gCenter = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zoneCenter, 0)], &geom, &elem);
  for (int mu=0; mu<NDIM; mu++)
  {
    uConCenter[mu] = elem.uCon[mu];
  }

  /* uConLeft */
  setGridZone(iTile, jTile,
              iInTile-1, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zone);
  getXCoords(&zone, CENTER, XCoords);
  setGeometry(XCoords, &geom); gLeft = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
  for (int mu=0; mu<NDIM; mu++)
  {
    uConLeft[mu] = elem.uCon[mu];
  }

  /* uConRight */
  setGridZone(iTile, jTile,
              iInTile+1, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zone);
  getXCoords(&zone, CENTER, XCoords);
  setGeometry(XCoords, &geom); gRight = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
  for (int mu=0; mu<NDIM; mu++)
  {
    uConRight[mu] = elem.uCon[mu];
  }

  for (int mu=0; mu<NDIM; mu++)
  {
    graduConVis[mu] = 
      slopeLimitedDerivative(uConLeft[mu],
                             uConCenter[mu],
                             uConRight[mu])/zoneCenter.dX1;
  }

  graduConHigherOrderTerm1Vis[0] = 
    slopeLimitedDerivative(gLeft*uConLeft[1],
                           gCenter*uConCenter[1],
                           gRight*uConRight[1])/zoneCenter.dX1;

  #if (COMPUTE_DIM==2)
    /* uConLeft */
    setGridZone(iTile, jTile,
                iInTile, jInTile-1,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);
    getXCoords(&zone, CENTER, XCoords);
    setGeometry(XCoords, &geom); gLeft = sqrt(-geom.gDet);
    setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
    for (int mu=0; mu<NDIM; mu++)
    {
      uConLeft[mu] = elem.uCon[mu];
    }

    /* uConRight */
    setGridZone(iTile, jTile,
                iInTile, jInTile+1,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);
    getXCoords(&zone, CENTER, XCoords);
    setGeometry(XCoords, &geom); gRight = sqrt(-geom.gDet);
    setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
    for (int mu=0; mu<NDIM; mu++)
    {
      uConRight[mu] = elem.uCon[mu];
    }

    for (int mu=0; mu<NDIM; mu++)
    {
      graduConVis[mu+NDIM] = 
        slopeLimitedDerivative(uConLeft[mu],
                               uConCenter[mu],
                               uConRight[mu])/zoneCenter.dX2;
    }

    graduConHigherOrderTerm1Vis[1] = 
      slopeLimitedDerivative(gLeft*uConLeft[2],
                             gCenter*uConCenter[2],
                             gRight*uConRight[2])/zoneCenter.dX2;

  #endif

#endif /* CONDUCTION */
}
