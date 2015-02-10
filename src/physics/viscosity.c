#include "physics.h"

void addViscositySourceTermsToResidual
(
  const REAL primTile[ARRAY_ARGS TILE_SIZE],
  ARRAY(primGlobal), ARRAY(primHalfStepGlobal), ARRAY(primOldGlobal),
  ARRAY(connectionGlobal),
  ARRAY(graduConVisGlobal), 
  ARRAY(graduConHigherOrderTerm1VisGlobal),
  ARRAY(graduConHigherOrderTerm2VisGlobal),
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
      REAL graduConHigherOrderTerm2Vis[COMPUTE_DIM];
      computeViscositySpatialGradientTerms
        (primTile, 
         iTile, jTile, 
         iInTile, jInTile,
         X1Start, X2Start,
         X1Size, X2Size, 
         graduConVis, 
         graduConHigherOrderTerm1Vis,
         graduConHigherOrderTerm2Vis
        );
      
      for (int mu=0; mu<NDIM; mu++)
      {
        INDEX_PETSC(graduConVisGlobal, &zoneCenter, mu) = graduConVis[mu];
      }

      INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 0) = 
        graduConHigherOrderTerm1Vis[0];

      INDEX_PETSC(graduConHigherOrderTerm2VisGlobal, &zoneCenter, 0) = 
        graduConHigherOrderTerm2Vis[0];

      #if (COMPUTE_DIM==2)
        for (int mu=0; mu<NDIM; mu++)
        {
          INDEX_PETSC(graduConVisGlobal, &zoneCenter, mu+NDIM) = graduConVis[mu+NDIM];
        }

        INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 1) = 
          graduConHigherOrderTerm1Vis[1];

        INDEX_PETSC(graduConHigherOrderTerm2VisGlobal, &zoneCenter, 1) = 
          graduConHigherOrderTerm2Vis[1];
      #endif
    } 
    else
    {
      #if (TIME_STEPPING==IMPLICIT)
        REAL graduConVis[COMPUTE_DIM*NDIM];
        REAL graduConHigherOrderTerm1Vis[COMPUTE_DIM];
        REAL graduConHigherOrderTerm2Vis[COMPUTE_DIM];
        computeViscositySpatialGradientTerms
          (primTile, 
           iTile, jTile, 
           iInTile, jInTile,
           X1Start, X2Start,
           X1Size, X2Size, 
           graduConVis,
           graduConHigherOrderTerm1Vis,
           graduConHigherOrderTerm2Vis
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
      REAL T    =   (ADIABATIC_INDEX-1.)
                  * INDEX_PETSC(primGlobal, &zoneCenter, UU)
                  / INDEX_PETSC(primGlobal, &zoneCenter, RHO);
      REAL TOld =   (ADIABATIC_INDEX-1.)
                  * INDEX_PETSC(primOldGlobal, &zoneCenter, UU)
                  / INDEX_PETSC(primOldGlobal, &zoneCenter, RHO);
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)
        REAL TCenter =   (ADIABATIC_INDEX-1.)
                       * elemCenter.primVars[UU]
                       / elemCenter.primVars[RHO];
      #elif (TIME_STEPPING == IMPLICIT)
        REAL TCenter = (T + TOld)/2.;
      #endif


      /* Higher order term 1 
       * (Psi+1/2/beta)*(sqrt(g)*u^\mu)_{,\mu} *
       * where beta = tauVis/2/eta */
      REAL beta       = elem.tauVis/(elem.eta*2.);
      REAL betaOld    = elemOld.tauVis/(elemOld.eta*2.);
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)

        REAL betaCenter = elemCenter.tauVis/(elemCenter.eta*2.);

      #elif (TIME_STEPPING == IMPLICIT)

        REAL betaCenter = 0.5*(  (elem.tauVis/(elem.eta*2.) )
                               + (elemOld.tauVis/(elemOld.eta*2.) )
                              );  
      #endif

      REAL dHigherOrderTerm1[NDIM];
      
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)

        dHigherOrderTerm1[0] = 
          (elemCenter.primVars[PSI]+0.5/betaCenter) 
        * sqrt(-geomCenter.gDet) 
        * (elem.uCon[0] - elemOld.uCon[0])/dt;
      
        dHigherOrderTerm1[1] = 
          (elemCenter.primVars[PSI]+0.5/betaCenter)
        * INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 0);

        dHigherOrderTerm1[2] = 0.; dHigherOrderTerm1[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm1[2] = 
            (elemCenter.primVars[PSI]+0.5/betaCenter)
          * INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 1);
        #endif

      #elif (TIME_STEPPING == IMPLICIT)

        dHigherOrderTerm1[0] = 
          0.5 * (elem.primVars[PSI] + elemOld.primVars[PSI]
		 +0.5/beta + 0.5/betaOld)
              * sqrt(-geomCenter.gDet) 
              * (elem.uCon[0] - elemOld.uCon[0])/dt;
      
        dHigherOrderTerm1[1] = 
          0.5*(  (  (elemOld.primVars[PSI] + 0.5/betaOld)
                  * INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 0)
                 )
               + 
                 (  (elem.primVars[PSI] + 0.5/beta)
                  * graduConHigherOrderTerm1Vis[0]
                 )
              );

        dHigherOrderTerm1[2] = 0.; dHigherOrderTerm1[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm1[2] = 
	    0.5*(  (  (elemOld.primVars[PSI] + 0.5/betaOld)
                  * INDEX_PETSC(graduConHigherOrderTerm1VisGlobal, &zoneCenter, 1)
                 )
               + 
                (  (elem.primVars[PSI] + 0.5/beta)
                  * graduConHigherOrderTerm1Vis[1]
                 )
	      );
        #endif

      #endif /* TIME_STEPPING options for higherOrderTerm1 */

      /* Higher order term 2 
       * (Psi*T)/(2*beta) * (beta*g*u^\mu/T)_{,\mu}*/

      REAL dHigherOrderTerm2[NDIM];
      
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)

        dHigherOrderTerm2[0] = 
          (elemCenter.primVars[PSI] * TCenter)
        / (2.*betaCenter) 
        * sqrt(-geomCenter.gDet) 
        * ((beta*elem.uCon[0]/T) - (betaOld*elemOld.uCon[0]/TOld))/dt;
      
        dHigherOrderTerm2[1] = 
          (elemCenter.primVars[PSI] * TCenter)
        / (2.*betaCenter) 
        * INDEX_PETSC(graduConHigherOrderTerm2VisGlobal, &zoneCenter, 0);

        dHigherOrderTerm2[2] = 0.; dHigherOrderTerm2[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm2[2] = 
            (elemCenter.primVars[PSI] * TCenter)
          / (2.*betaCenter) 
          * INDEX_PETSC(graduConHigherOrderTerm2VisGlobal, &zoneCenter, 1);
        #endif

      #elif (TIME_STEPPING == IMPLICIT)

        dHigherOrderTerm2[0] = 
          0.5*(  (elem.primVars[PSI] * T/(2.*beta) )
               + (elemOld.primVars[PSI] * TOld/(2.*betaOld) )
              )
        * sqrt(-geomCenter.gDet) 
        * ((beta*elem.uCon[0]/T) - (betaOld*elemOld.uCon[0]/TOld))/dt;


        dHigherOrderTerm2[1] = 
          0.5*(  (  elemOld.primVars[PSI] * TOld/(2.*betaOld) 
                  * INDEX_PETSC(graduConHigherOrderTerm2VisGlobal, &zoneCenter, 0)
                 )
               + 
                 (  elem.primVars[PSI] * T/(2.*beta)
                  * graduConHigherOrderTerm2Vis[0]
                 )
              );

        dHigherOrderTerm2[2] = 0.; dHigherOrderTerm2[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm2[2] = 
            0.5*(  (  elemOld.primVars[PSI] * TOld/(2.*betaOld) 
                    * INDEX_PETSC(graduConHigherOrderTerm2VisGlobal, &zoneCenter, 1)
                   )
                + 
                   (  elem.primVars[PSI] * T/(2.*beta)
                    * graduConHigherOrderTerm2Vis[1]
                   )
                );
        #endif

      #endif

      REAL higherOrderTerm1 = 0.;
      REAL higherOrderTerm2 = 0.;
      for (int mu=0; mu<NDIM; mu++)
      {
        higherOrderTerm1 += dHigherOrderTerm1[mu];
        higherOrderTerm2 += dHigherOrderTerm2[mu];
      }


      // Term g*(3*eta*b^mu b_\nu/b^2 u^\nu_{;\nu}+Psi)/tau
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
	      * (elem.uCon[mu] - elemOld.uCon[mu])/dt;
	    TargetPsi+=1.5*(bCov[alpha]*elem.bCon[1]*elem.eta*graduConVis[alpha]/bSqr
			    +bCovOld[alpha]*elemOld.bCon[1]*elemOld.eta/bSqrOld
			    *INDEX_PETSC(graduConVisGlobal, &zoneCenter, alpha));
            #if (COMPUTE_DIM==2)
              TargetPsi+=1.5*(bCov[alpha]*elem.bCon[2]*elem.eta*graduConVis[alpha]/bSqr
			      +bCovOld[alpha]*elemOld.bCon[2]*elemOld.eta/bSqrOld
			      *INDEX_PETSC(graduConVisGlobal, &zoneCenter, alpha));
            #endif
	  }
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
	for (int alpha=0; alpha<NDIM; alpha++)
          {
            for (int beta=0; beta<NDIM; beta++)
	      {
		for(int mu=0; mu<NDIM; mu++)
		  {
		    TargetPsi+=
		      1.5*INDEX_PETSC(connectionGlobal, &zoneCenter, GAMMA_UP_DOWN_DOWN(mu, alpha, beta))
		      *(elem.bCon[alpha]*bCov[mu]*elem.uCon[beta]/bSqr*elem.eta
			+elemOld.bCon[alpha]*bCovOld[mu]*elemOld.uCon[beta]/bSqrOld*elemOld.eta);
		  }
	      }
	  }
      #elif (TIME_STEPPING == EXPLICIT)
	REAL bSqr, bCov[NDIM];
        conToCov(elemCenter.bCon, &geomCenter, bCov);
        bSqr = getbSqr(&elemCenter, &geomCenter);
        for (int alpha=0; alpha<NDIM; alpha++)
          {
            for (int beta=0; beta<NDIM; beta++)
	      {
		for(int mu=0; mu<NDIM; mu++)
                  {
                    TargetPsi+=
                      3.*INDEX_PETSC(connectionGlobal, &zoneCenter, GAMMA_UP_DOWN_DOWN(mu, alpha, beta))
                      *(elemCenter.bCon[alpha]*bCov[mu]*elemCenter.uCon[beta]/bSqr*elemCenter.eta);
                  }
	      }
	  }
      #endif	
      }

      #if (TIME_STEPPING == EXPLICIT)

        INDEX_PETSC(residualGlobal, &zoneCenter, PSI) += 
          (- higherOrderTerm1 + higherOrderTerm2
           + g*( elemCenter.primVars[PSI]
                 + TargetPsi
               )/elemCenter.tauVis
          )/norm;

      #elif (TIME_STEPPING == IMEX || TIME_STEPPING == IMPLICIT)
  
        INDEX_PETSC(residualGlobal, &zoneCenter, PSI) += 
          (- higherOrderTerm1 + higherOrderTerm2
           + g*( 0.5*(elem.primVars[PSI] + elemOld.primVars[PSI])
                + TargetPsi
               )/elemCenter.tauVis
          )/norm;

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
  REAL graduConHigherOrderTerm1Vis[COMPUTE_DIM],
  REAL graduConHigherOrderTerm2Vis[COMPUTE_DIM]
)
{
#if (VISCOSITY)
  REAL XCoords[NDIM];
  struct gridZone zone, zoneCenter;
  struct geometry geom;
  struct fluidElement elem;
  REAL uConLeft[NDIM], uConCenter[NDIM], uConRight[NDIM];
  REAL gLeft, gCenter, gRight;
  REAL betaLeft, betaCenter, betaRight; /* Israel-Stewart's Beta */

  setGridZone(iTile, jTile,
              iInTile, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zoneCenter);

  REAL TCenter  =   (ADIABATIC_INDEX-1.)
                  * primTile[INDEX_TILE(&zoneCenter, UU)]
                  / primTile[INDEX_TILE(&zoneCenter, RHO)];
  REAL TLeftX1  =   (ADIABATIC_INDEX-1.)
                  * primTile[INDEX_TILE_MINUS_ONE_X1(&zoneCenter, UU)]
                  / primTile[INDEX_TILE_MINUS_ONE_X1(&zoneCenter, RHO)];
  REAL TRightX1 =  (ADIABATIC_INDEX-1.)
                  * primTile[INDEX_TILE_PLUS_ONE_X1(&zoneCenter, UU)]
                  / primTile[INDEX_TILE_PLUS_ONE_X1(&zoneCenter, RHO)];
  #if (COMPUTE_DIM==2)
    REAL TLeftX2  =   (ADIABATIC_INDEX-1.)
                    * primTile[INDEX_TILE_MINUS_ONE_X2(&zoneCenter, UU)]
                    / primTile[INDEX_TILE_MINUS_ONE_X2(&zoneCenter, RHO)];
    REAL TRightX2 =   (ADIABATIC_INDEX-1.)
                    * primTile[INDEX_TILE_PLUS_ONE_X2(&zoneCenter, UU)]
                    / primTile[INDEX_TILE_PLUS_ONE_X2(&zoneCenter, RHO)];
  #endif


  /* uConCenter */
  getXCoords(&zoneCenter, CENTER, XCoords);
  setGeometry(XCoords, &geom); gCenter = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zoneCenter, 0)], &geom, &elem);
  betaCenter = elem.tauVis/(elem.eta*2.);
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
  betaLeft = elem.tauVis/(elem.eta*2.);
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
  betaRight = elem.tauVis/(elem.eta*2.);
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

  graduConHigherOrderTerm2Vis[0] = 
    slopeLimitedDerivative
      (betaLeft*gLeft*uConLeft[1]/TLeftX1,
       betaCenter*gCenter*uConCenter[1]/TCenter,
       betaRight*gRight*uConRight[1]/TRightX1
      )/zoneCenter.dX1;

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
    betaLeft = elem.tauVis/(elem.eta*2.);
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
    betaRight = elem.tauVis/(elem.eta*2.);
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

    graduConHigherOrderTerm2Vis[1] = 
      slopeLimitedDerivative
        (betaLeft*gLeft*uConLeft[2]/TLeftX2,
         betaCenter*gCenter*uConCenter[2]/TCenter,
         betaRight*gRight*uConRight[2]/TRightX2
        )/zoneCenter.dX2;

  #endif

#endif /* CONDUCTION */
}
